#!/usr/bin/env python3

import argparse
import subprocess
from collections import defaultdict
import pysam


parser = argparse.ArgumentParser(
    description="4D heterozygosity (π proxy) for single diploid sample"
)

parser.add_argument("-r", "--reference", required=True)
parser.add_argument("-g", "--gff", required=True)
parser.add_argument("-v", "--vcf", required=True)

args = parser.parse_args()

genome_fa = args.reference
gff_file = args.gff
vcf_file = args.vcf


fourfold = {"GC", "GG", "CC", "AC", "GT", "CG", "CT", "TC"}


def fetch_seq(chrom, start, end):
    cmd = f"samtools faidx {genome_fa} {chrom}:{start+1}-{end}"
    out = subprocess.check_output(cmd, shell=True).decode().splitlines()
    return "".join(out[1:]).upper()


genes = defaultdict(list)

with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue

        chrom, _, feat, start, end, _, strand, phase, attrs = line.strip().split("\t")

        if feat != "CDS":
            continue

        raw_parent = next(
            (x.split("=")[1] for x in attrs.split(";") if x.startswith("Parent=")),
            None
        )

        # remove chrom prefix (e.g. chrom1-gene123 → gene123)
        parent = raw_parent.split("-", 1)[-1] if raw_parent else None

        genes[parent].append(
            (chrom, int(start) - 1, int(end), strand, int(phase))
        )


fourD_sites = []

for blocks in genes.values():

    blocks.sort(key=lambda x: x[1])

    strand = blocks[0][3]

    if strand == "-":
        blocks = blocks[::-1]

    seq = []
    coord = []

    for chrom, start, end, strand, phase in blocks:

        frag = fetch_seq(chrom, start, end)

        if strand == "-":
            frag = frag[::-1].translate(str.maketrans("ACGT", "TGCA"))

        frag = frag[phase:]

        for i, base in enumerate(frag):
            seq.append(base)
            coord.append((chrom, start + i))

    seq = seq[:len(seq) - len(seq) % 3]

    for i in range(0, len(seq), 3):
        codon = "".join(seq[i:i+3])

        if "N" in codon:
            continue

        if codon[:2] in fourfold:
            fourD_sites.append(coord[i + 2])


vcf = pysam.VariantFile(vcf_file)

fourD_set = set(fourD_sites)

hets = 0

for record in vcf.fetch():

    key = (record.chrom, record.pos - 1)

    if key not in fourD_set:
        continue

    gt = next(iter(record.samples.values()))["GT"]

    if gt is not None and len(gt) == 2 and gt[0] != gt[1]:
        hets += 1


pi = hets / len(fourD_sites) if fourD_sites else 0

print("=== 4D Heterozygosity (π proxy) ===")
print(f"Reference:     {genome_fa}")
print(f"GFF:           {gff_file}")
print(f"VCF:           {vcf_file}")
print(f"4D sites:      {len(fourD_sites)}")
print(f"Heterozygotes: {hets}")
print(f"π (proxy):     {pi:.6e}")
