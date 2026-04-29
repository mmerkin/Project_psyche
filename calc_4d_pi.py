#!/usr/bin/env python3

import argparse
from collections import defaultdict
import pysam
import time

# Parse inputs

parser = argparse.ArgumentParser(
    description="CDS extraction + flexible 0/2/3/4-fold VCF output"
)

parser.add_argument("-r", "--reference", required=True)
parser.add_argument("-g", "--gff", required=True)
parser.add_argument("-v", "--vcf", required=True)
parser.add_argument("-o", "--out_prefix", required=True)
parser.add_argument("-f", "--fold", required=True,
                    help="Comma-separated folds to output, e.g. 0,2,4")

args = parser.parse_args()

genome_fa = args.reference
gff_file = args.gff
vcf_file = args.vcf
prefix = args.out_prefix

requested_folds = set(args.fold.split(","))


# Set output files

cds_file = f"{prefix}_cds.fa"
cds_out = open(cds_file, "w")

vcf = pysam.VariantFile(vcf_file)

fasta = pysam.FastaFile(genome_fa)

RC = str.maketrans("ACGT", "TGCA")


# Degeneracy table taken from Mackintosh 2019

degen_dict = {
    'ttt': '002', 'ttc': '002', 'tta': '202', 'ttg': '202',
    'tct': '004', 'tcc': '004', 'tca': '004', 'tcg': '004',
    'tat': '002', 'tac': '002', 'taa': '022', 'tag': '002',
    'tgt': '002', 'tgc': '002', 'tga': '020', 'tgg': '000',
    'ctt': '004', 'ctc': '004', 'cta': '204', 'ctg': '204',
    'cct': '004', 'ccc': '004', 'cca': '004', 'ccg': '004',
    'cat': '002', 'cac': '002', 'caa': '002', 'cag': '002',
    'cgt': '004', 'cgc': '004', 'cga': '204', 'cgg': '204',
    'att': '003', 'atc': '003', 'ata': '003', 'atg': '000',
    'act': '004', 'acc': '004', 'aca': '004', 'acg': '004',
    'aat': '002', 'aac': '002', 'aaa': '002', 'aag': '002',
    'agt': '002', 'agc': '002', 'aga': '202', 'agg': '202',
    'gtt': '004', 'gtc': '004', 'gta': '004', 'gtg': '004',
    'gct': '004', 'gcc': '004', 'gca': '004', 'gcg': '004',
    'gat': '002', 'gac': '002', 'gaa': '002', 'gag': '002',
    'ggt': '004', 'ggc': '004', 'gga': '004', 'ggg': '004'
}


genes = defaultdict(list)


start_time = time.time()


# Convert the gff file to a bed with just CDS locations 

with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split()
        if len(parts) < 9: # Skips malformed lines during testing, probably unnecessary now
            continue

        chrom, _, feat, start, end, _, strand, phase, attrs = parts[:9]

        if feat.upper() != "CDS": # Skips gene, mRNA etc lines
            continue

        parent = next(
            (x.split("=", 1)[1] for x in attrs.split(";") if x.startswith("Parent=")), # Extract the gene name from the attrs sections (after parent=)
            None
        )

        if not parent:
            continue

        genes[parent].append(
            (chrom, int(start) - 1, int(end), strand, int(phase)) # Converts 1-based GFF to 0-based bed
        )


# Create writers for output vcf files

writers = {}

for f in requested_folds:
    writers[f] = pysam.VariantFile(
        f"{prefix}.{f}d.vcf", "w", header=vcf.header
    )


# Keep track of all sites and heterozygous sites

sites = {
    "0": set(),
    "2": set(),
    "3": set(),
    "4": set()
}

heterozygous_sites = {
    "0": 0,
    "2": 0,
    "3": 0,
    "4": 0
}



# Find CDS, split it into codons, compare the third position to a degen table and add to sites

for parent, blocks in genes.items():

    blocks.sort(key=lambda x: (x[0], x[1]))

    seq = []
    coord = []

    for chrom, start, end, strand, phase in blocks:

        frag = fasta.fetch(chrom, start, end).upper()

        if strand == "-":
            frag = frag[::-1].translate(RC) # Reverse compliment genes on minus strand

        for i, base in enumerate(frag):
            seq.append(base)
            coord.append((chrom, start + i))


    gene_seq = "".join(seq)

    if gene_seq:
        cds_out.write(f">{parent}\n") 
        for i in range(0, len(gene_seq), 60):     # write CDS fasta with a wrap of 60 bases per line
            cds_out.write(gene_seq[i:i+60] + "\n")

    seq = seq[:len(seq) - (len(seq) % 3)] # Remove bases left if they are in an incomplete codon
    coord = coord[:len(coord) - (len(coord) % 3)] # Same for coordinates

    for i in range(0, len(seq), 3):

        codon = "".join(seq[i:i+3]).lower()

        if "n" in codon:
            continue

        site = coord[i + 2] # Get the coordinate of the third position of the codon

        deg = degen_dict.get(codon)
        if deg is None:
            continue

        third = deg[2] # Use the look-up table to see what the degeneracy of the third position is

        if third in sites:
            sites[third].add(site)


# Make VCF files

written = {f: 0 for f in requested_folds}

for record in vcf.fetch():

    key = (record.chrom, record.pos - 1) # Convert 1-based vcf variants to 0-based

    for f in requested_folds: # Add variants to a different vcf depending on fold
        if key in sites[f]:
            writers[f].write(record)
            written[f] += 1
            
            for sample in record.samples.values():
                if sample["GT"] is not None and len(set(sample["GT"])) > 1:  # Heterozygous (GT not equal to a single value)
                    heterozygous_sites[f] += 1
                    break


for w in writers.values():
    w.close()
cds_out.close()


end_time = time.time() 
elapsed_time = end_time - start_time


# Summary to be printed to terminal

print("Complete")

for f in sorted(requested_folds):
    print(f"{f}D sites:", len(sites[f]))
    print(f"{f}D VCF records:", written[f])
    heterozygosity = heterozygous_sites[f] / written[f] if written[f] > 0 else 0
    print(f"{f}D Heterozygosity: {heterozygosity:.4f}")

print(f"Total execution time: {elapsed_time:.2f} seconds")
