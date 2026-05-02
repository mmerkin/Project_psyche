# Calling variants

Code for calling variants

```bash
run_sample() {
  reads="$1"
  threads="$2"
  
  
  j=$((threads / 2))
  J=$((threads - j))

  species=$(basename "$reads" | cut -d'_' -f1,2)
  genome=$(ls Genomes/${species}.GCA_*/*.fa)
  pbmm2 index "$genome" "${genome%.*}.mmi"
  pbmm2 align "${genome%.*}.mmi" "$reads" "${species}.sorted.bam" --rg "@RG\tID:${species}_pb_1\tSM:${species}\tPL:PACBIO\tLB:${species}\tPU:${species}_run1" --sort -j "$j" -J "$J"
  samtools depth -a -Q 20 -q 30 "${species}.sorted.bam" > "${species}.depth.txt"
  mean_cov=$(awk '{sum+=$3} END {print sum/NR}' "${species}.depth.txt")
  max_cov=$(awk -v m="$mean_cov" 'BEGIN {print 2*m}')

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
}

export -f run_sample


parallel -j 4 run_sample {} 12 ::: HiFi_reads/*.fastq.gz

```


## Genomescope

```bash
samtools view -bS -L Aricia_artaxerxes_autosomes.bed Aricia_artaxerxes.sorted.bam > Aricia_artaxerxes_autosomes.bam
samtools bam2fq Aricia_artaxerxes_autosomes.bam > Aricia_artaxerxes_autosomes.fastq; bgzip Aricia_artaxerxes_autosomes.fastq


FastK -v -t1 -k31 ../Aricia_artaxerxes_HiFi.fastq.gz -NTable
Histex -G Table | ~/apps/GENESCOPE.FK/GeneScopeFK.R -o Output -k 31







## Longshot

Old code for a variant caller that is no longer being used 

```bash

#mkdir -p logs

cut -f1 "$genome.fai" | parallel -j 8 '
longshot --bam '"${species}"'.sorted.bam --ref '"$genome"' \
--out out_{}.vcf --region {} \
> logs/{}.log 2>&1'

parallel -j 8 bgzip {} ::: out_*.vcf
parallel -j 8 tabix {} ::: out_*.vcf.gz

bcftools concat -a out_*.vcf.gz -Oz -o merged.vcf.gz
bcftools index merged.vcf.gz
```
