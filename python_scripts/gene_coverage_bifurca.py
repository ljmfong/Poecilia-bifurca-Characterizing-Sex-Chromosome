'''
This script calculates Male to female read depth ratio 
Input: male.cds.ldepth female.cds.ldepth
Output: M_F.cds.coverage.ratio.csv
Usage: python gene_coverage.py gene_boundary.bed male.cds.ldepth female.cds.ldepth
Date: Jan 13, 2021
'''
import sys,os
import pandas as pd
import math

def get_gene_boundary(bed):
	with open(bed,'r') as inf:
		header = inf.readline()
		for line in inf:
			line = line.rstrip("\n")
			chrom, start, end, gene_IDs = line.split("\t") 
			#chrom, start, end, gene_IDs = line.split("\t") #Depends on how your bed file is separated
			yield chrom, start, end, gene_IDs


def get_snps(file,chrom,start,end, gene_IDs):
	# CHROM	POS	male_sum	SUMSQ_DEPTH	female_sum	SUMSQ_DEPTH.1	MF_ratio
	df = pd.read_csv(file, sep = "\t")
	subset = df[(df["CHROM"]== str(chrom)) & (df["POS"] <= int(end)) & (df["POS"] > int(start))]
	if len(subset) == 0:
		return None 
	else:
		return subset

output_file = sys.argv[3]

with open(output_file, 'w') as outf:
    outf.write('\t'.join(["CHROM", "START", "END", "Gene_ID", "POS", "MF_ratio", "FST"]) + "\n")

    for chrom, start, end, gene_IDs in get_gene_boundary(sys.argv[1]):
        matching_snps = get_snps(sys.argv[2], chrom, start, end, gene_IDs)
        if matching_snps is not None:
            for index, row in matching_snps.iterrows():
                outf.write('\t'.join([chrom, str(start), str(end), gene_IDs, str(row["POS"]), str(row["MF_ratio"])]) + "\n")
