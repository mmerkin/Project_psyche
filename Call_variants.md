# Calling variants

Code for calling variants

```bash
run_sample() {
  reads="$1"
  output="$2"
  threads="$3"
  
  
  j=$((threads / 2))
  J=$((threads - j))

  species=$(basename "$reads" | cut -d'_' -f1,2)
  genome=$(ls Genomes/${species}.GCA_*/*.fa)
  pbmm2 index "$genome" "${genome%.*}.mmi"
  pbmm2 align "${genome%.*}.mmi" "$reads" "$output/${species}.sorted.bam" --rg "@RG\tID:${species}_pb_1\tSM:${species}\tPL:PACBIO\tLB:${species}\tPU:${species}_run1" --sort -j "$j" -J "$J"
  samtools depth -a -Q 20 -q 30 "$output/${species}.sorted.bam" > "$output/${species}.depth.txt"
}

export -f run_sample


parallel -j 4 run_sample {} 12 ::: HiFi_reads/*.fastq.gz

```

```bash

  bcftools mpileup -f "$genome" \
    -a FORMAT/QS,FORMAT/AD,FORMAT/DP,INFO/AD \
    -B --min-MQ 30 --min-BQ 20 \
    "${species}.sorted.bam" | \
  bcftools call -m --threads "$threads" --ploidy 2 -a GQ,GP \
    -Oz -o "${species}.bcftools.vcf.gz"

  bcftools view -v snps -m2 -M2 \
    "${species}.bcftools.vcf.gz" \
    -Oz -o "${species}.bcftools.snps.vcf.gz"
    
  bcftools filter \
  -i "FORMAT/DP>=6 && FORMAT/DP<=${max_cov} && FORMAT/GQ>20" \
  "${species}.bcftools.snps.vcf.gz" -Oz -o "${species}.bcftools.snps.filtered.vcf.gz"

  tabix "${species}.bcftools.snps.filtered.vcf.gz"


```


## Genomescope

```bash
samtools view -bS -L Aricia_artaxerxes_autosomes.bed Aricia_artaxerxes.sorted.bam > Aricia_artaxerxes_autosomes.bam
samtools bam2fq Aricia_artaxerxes_autosomes.bam > Aricia_artaxerxes_autosomes.fastq; bgzip Aricia_artaxerxes_autosomes.fastq


FastK -v -t1 -k31 ../Aricia_artaxerxes_HiFi.fastq.gz -NTable
Histex -G Table | ~/apps/GENESCOPE.FK/GeneScopeFK.R -o Output -k 31

```

## Python

```bash
cat Heterozygosity/Agriades_optilete_autosomes.heterozygosity.txt | awk '$4 >= 6000 {het += $6; cov += $4} END {print het/cov}'
```



## Longshot



```bash

threads=40
for bam in bams/*.bam; do
species=$(basename "$bam" | cut -d'.' -f1)
genome=$(ls Genomes/${species}.GCA_*/*.fa)

mkdir -p "Variants/$species/logs"

cut -f1 "$genome.fai" | parallel -j "$threads" \
longshot --bam "$bam" --ref "$genome" \
--out Variants/"$species"/"${species}_chr_{}".vcf --region {} \
">" Variants/"$species"/logs/"${species}_chr_{}".log "2>&1"

parallel -j "$threads" bgzip ::: Variants/"$species"/${species}_chr_*.vcf
parallel -j "$threads" tabix ::: Variants/"$species"/${species}_chr_*.vcf.gz

bcftools concat -a Variants/"$species"/${species}_chr_*.vcf.gz -Oz -o Variants/"$species"/${species}_merged.vcf.gz
bcftools index Variants/"$species"/${species}_merged.vcf.gz
done
```

```bash

while read species; do
bam=$species.sorted.bam
samtools depth bams/$bam > depths/$bam.depth.txt
awk '{sum+=$3} END { print "Average = ",sum/NR}' depths/$bam.depth.txt > depths/$bam.cov
done < Species_list.txt



cat Species_list.txt | parallel 'samtools depth bams/{}.sorted.bam | awk '\''{sum+=$3} END {print "Average = ",sum/NR}'\'' > depths/{}.sorted.bam.cov'

cat Species_list.txt | parallel scripts/updated_filter_variants.py --vcf Variants/{}/{}_merged.vcf.gz --cov depths/{}.sorted.bam.cov --o Filtered_variants/{}

for i in Filtered_variants/*/*vcf.gz; do tabix $i; done

cat Species_list.txt | parallel scripts/updated_calculate_genome_wide_heterozygosity.py --vcf Filtered_variants/{}/{}_merged.filtered.vcf.gz --bam bams/{}.sorted.bam --cov depths/{}.sorted.bam.cov --o Heterozygosity/{}

## Needs to be saved to an output file and also fixed to account for uncallable sites

#cat Species_list.txt | parallel scripts/calc_fold_pi.py -r $(find Genomes/{}* -name "*.fa") -g $(find Annotations/{}* -name "*.gff") -v Variants/{}/{}_merged.vcf.gz -o PI/{} -f 0,4



## Issues to fix:
### Running the same command again appends rather than creating a new file (extract chromosomes for het)
### Some individuals have 0 ROHs (especially Z), which breaks script updated_ROH_size_categories.py
### Should also calculate 4d pi 

species_list=($(cat Species_list.txt))
z_list=($(cat Z_list.txt))


num_species=${#species_list[@]}
for ((i=0; i<$num_species; i++)); do
species="${species_list[$i]}"
Z="${z_list[$i]}"


bash scripts/extract_chromosomes_for_het.sh $(find Genomes/$species*/ -name "*autosomes.txt") "Heterozygosity/$species/${species}_merged.filtered.heterozygosity.txt" "Heterozygosity/$species/${species}_autosomes.heterozygosity.txt"
bash scripts/extract_chromosomes_for_het.sh "$Z" "Heterozygosity/$species/${species}_merged.filtered.heterozygosity.txt" "Heterozygosity/$species/${species}_Z2.heterozygosity.txt"

python3 scripts/runs_of_homozygosity.py --het Heterozygosity/$species/${species}_autosomes.heterozygosity.txt --o ROH/$species
python3 scripts/runs_of_homozygosity.py --het Heterozygosity/$species/${species}_Z2.heterozygosity.txt --o ROH/$species

python3 scripts/updated_ROH_size_categories.py --het ROH/$species/${species}_autosomes.heterozygosity.insideROH.txt --o ROH/$species
python3 scripts/updated_ROH_size_categories.py --het ROH/$species/${species}_Z2.heterozygosity.insideROH.txt --o ROH/$species

cat Heterozygosity/$species/${species}_autosomes.heterozygosity.txt | awk '$4 >= 6000 {het += $6; cov += $4} END {print het/cov}'
cat Heterozygosity/$species/${species}_Z2.heterozygosity.txt | awk '$4 >= 6000 {het += $6; cov += $4} END {print het/cov}'

done

```
## Old code for a variant caller that is no longer being used 

```bash
bcftools view -v snps merged.vcf.gz -Ou | \
bcftools filter -i "QUAL>=15 && FORMAT/GQ>=20" \
-Oz -o filtered.vcf.gz


vcftools --gzvcf filtered.vcf.gz --SNPdensity 10000 --out snp_density

awk -v maxdp="$MAX_DP" '{
  if ($3 >= 6 && $3 <= maxdp)
    print $1"\t"$2-1"\t"$2
}' depth.txt > callable.bed

sort -k1,1 -k2,2n callable.bed | bedtools merge > callable.merged.bed

bedtools makewindows -g "$genome.fai" -w 10000 > windows.bed

bedtools coverage -a windows.bed -b callable.merged.bed > callable_per_window.txt

paste snp_density.snpden callable_per_window.txt | \
awk '{
  chrom=$1
  start=$2
  end=$3
  snps=$4
  callable=$NF

  if (callable > 0) {
    norm = snps * (10000 / callable)
  } else {
    norm = "NA"
  }

  print chrom, start, end, snps, callable, norm
}' > normalized_snp_density.txt
```
