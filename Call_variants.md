# Calling variants

Code for calling variants

```bash
genome=GCA_937612035.1.fa
reads=Aricia_artaxerxes_HiFi.fastq.gz
species=Aricia_artaxerxes
threads=32

pbmm2 index "$genome" "${genome%.*}.mmi"
pbmm2 align "${genome%.*}.mmi" "$reads" "${species}.sorted.bam" --rg "@RG\tID:${species}_pb_1\tSM:${species}\tPL:PACBIO\tLB:${species}\tPU:${species}_run1" --sort -j 32 -J 16

bcftools mpileup -f $genome -a "FORMAT/QS,FORMAT/AD,FORMAT/DP,INFO/AD" -B --min-MQ 30 --min-BQ 20 ${species}.sorted.bam | \
bcftools call -m --threads $threads --ploidy 2 -a GQ,GP -Oz -o ${base}.bcftools.vcf.gz

bcftools view -v snps -m2 -M2 "${species}.bcftools.vcf.gz" -Oz -o "${species}.bcftools.snps.vcf.gz" && \
tabix "${species}.biallelic.snps.vcf.gz"

```

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
