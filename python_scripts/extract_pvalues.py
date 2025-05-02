#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
'''
'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
from collections import Counter
import vcf
import time
from scipy import stats as stats

# Use R via rpy2
# import rpy2.robjects as R
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
					help="consistent file")
parser.add_argument("outfile", type=str,
					help="out unadjusted pvalues")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	dicty = {}
	new_dicty = {}

	with open(args.infile, "r") as infile:
		next(infile)
		for line in infile:
			line = line.rstrip()
			scaffold = line.split("\t")[0]
			position = line.split("\t")[1]
			format = line.split("\t")[2]
			samples = line.split("\t")[3:6]
			for sample in samples:
				genotype = sample.split(":")[0]
				if genotype == "0/1":
					ids = scaffold + "_" + position + "_" + format + "_" + sample
					pvalue = sample.split(":")[4]
					dicty[ids] = pvalue
	
	with open(args.outfile, "w") as out:
		for entry in dicty:
			out.write(entry)
			out.write("\t")
			out.write(str(dicty[entry]))
			out.write("\n")

# for t in AUTO_ASE_snps:
#     AUTO_ASE_snps[t] = fdr_snps(AUTO_ASE_snps[t])


if __name__ == '__main__':
	main()
