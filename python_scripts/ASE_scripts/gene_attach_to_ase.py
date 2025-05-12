'''
This script outputs the gene ID for your ASE files
Input: BED file; ASE files (comma delimited)
Output: 
Usage: python gene_attach_to_ase.py gene_boundary.bed ase_file output_file
Date: Jan 13, 2021
'''
import sys,os
import pandas as pd
import math

# Load BED file (fileA)
def get_gene_boundary(bed):
    with open(bed, 'r') as inf:
        header = inf.readline()  # skip header
        for line in inf:
            line = line.rstrip("\n")
            chrom, start, end, gene_id = line.split("\t")
            yield chrom, int(start), int(end), gene_id

# Load CSV file (fileB)
df_b = pd.read_csv(sys.argv[2])  # columns: count, scaffold, position, average_majallele_fraction, chromo

# Output file
with open(sys.argv[3], 'w') as outf:
    outf.write('\t'.join(["count", "scaffold", "position", "average_majallele_fraction", "chromo", "Geneid"]) + '\n')

    # Loop over BED rows
    for chrom, start, end, gene_id in get_gene_boundary(sys.argv[1]):
        subset = df_b[(df_b["scaffold"] == chrom) & 
                      (df_b["position"] >= start) & 
                      (df_b["position"] <= end)]
        
        for _, row in subset.iterrows():
            fields = [
                str(row["count"]),
                row["scaffold"],
                str(row["position"]),
                str(row["average_majallele_fraction"]),
                row["chromo"],
                gene_id
            ]
            outf.write('\t'.join(fields) + '\n')

