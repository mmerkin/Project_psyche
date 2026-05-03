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

Old code for a variant caller that is no longer being used 

```bash

mkdir -p logs

samtools depth -a -q 20 sample.bam > depth.txt
MEAN_DP=$(awk '{sum+=$3} END {print sum/NR}' depth.txt)
MAX_DP=$((MEAN_DP * 2))

cut -f1 "$genome.fai" | parallel -j 8 '
longshot --bam '"${species}"'.sorted.bam --ref '"$genome"' \
--min_cov 6 --max_cov '"$MAX_DP"' --min_mapq 20 \
--out out_{}.vcf --region {} \
> logs/{}.log 2>&1'

parallel -j 8 bgzip {} ::: out_*.vcf
parallel -j 8 tabix {} ::: out_*.vcf.gz

bcftools concat -a out_*.vcf.gz -Oz -o merged.vcf.gz
bcftools index merged.vcf.gz

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
