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
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
					help="A vcf file")
parser.add_argument("out_density", type=str,
					help="outfile with info about the ref_allele_fraction for R plotting")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def binom_test_ase(rd,ad):
	'''
    RD: Depth of reference-supporting bases (reads1)
    AD: Depth of variant-supporting bases (reads2)
    This is the chance of observing either more successes, or fewer successes, in all trials.
    '''
    # Two tailed binomial test
	return stats.binom_test(rd,(rd+ad),p=0.5)
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	count = 0
	with open(args.out_density, "w") as out_density:
		out_density.write("count,scaffold,position,average_majallele_fraction\n")
		with open(args.infile, "r") as vcf:
			for line in vcf:
				if line.startswith("#"):
					pass
				if line.startswith("LG"):
					count_hetero = 0
					sum_fraction = 0
					line = line.rstrip()
					scaffold = line.split()[0]
					position = line.split()[1]
					samples = line.split()[9:12]
					for sample in samples:
						genotype = sample.split(":")[0]
						rd_value = sample.split(":")[4]
						ad_value = sample.split(":")[5]
						if rd_value == '.' or ad_value == '.':
							continue
						
						RD = float(rd_value)
						AD = float(ad_value)

						if RD > AD:
							maj_allele_fraction = (float(RD))/((float(RD))+(float(AD)))
						else:
							maj_allele_fraction = (float(AD))/((float(RD))+(float(AD)))
						if genotype == "0/1":
							count_hetero += 1
							sum_fraction += maj_allele_fraction
							#print(sum_fraction)
					if count_hetero >= 1:
						count += 1
						out_density.write(str(count))
						out_density.write(",")
						out_density.write(scaffold)
						out_density.write(",")
						out_density.write(str(position))
						average_fraction = sum_fraction/count_hetero
						out_density.write(",")
						out_density.write(str(average_fraction))
						out_density.write("\n")

if __name__ == '__main__':
	main()
