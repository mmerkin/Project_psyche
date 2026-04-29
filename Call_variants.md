# Calling variants

Code for calling variants

```bash
genome=GCA_937612035.1.fa
reads=Aricia_artaxerxes_HiFi.fastq.gz
species=Aricia_artaxerxes

pbmm2 index "$genome" "${genome%.*}.mmi"
pbmm2 align "${genome%.*}.mmi" "$reads" "${species}.sorted.bam" --rg "@RG\tID:${species}_pb_1\tSM:${species}\tPL:PACBIO\tLB:${species}\tPU:${species}_run1" --sort -j 32 -J 16

mkdir -p logs

cut -f1 "$genome.fai" | parallel -j 8 '
longshot --bam '"${species}"'.sorted.bam --ref '"$genome"' \
--out out_{}.vcf --region {} \
> logs/{}.log 2>&1'

bcftools concat -a out_*.vcf -Oz -o merged.vcf.gz
bcftools index merged.vcf.gz


bcftools view -v snps -m2 -M2 "${species}.vcf" -Oz -o "${species}.biallelic.snps.vcf.gz" && \
bcftools index "${species}.biallelic.snps.vcf.gz"
```
