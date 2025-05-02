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
parser.add_argument("outfile", type=str,
					help="outfile with info about ase - fabian parameters")
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
		out_density.write("count,scaffold,position,fraction_s1,fraction_s2,fraction_s3\n")
		with open(args.outfile, "w") as out:
			out.write("scaffold\tposition\tformat\t")
			with open(args.infile, "r") as vcf:
				for line in vcf:
					if line.startswith("##"):
						pass
					if line.startswith("#CHROM"):
						sample1 = line.split()[9]
						sample2 = line.split()[10]
						sample3 = line.split()[11]
						out.write(sample1)
						out.write("\t")
						out.write(sample2)
						out.write("\t")
						out.write(sample3)
						out.write("\n")
					if line.startswith("Scaff"):
						count += 1
						line = line.rstrip()
						scaffold = line.split()[0]
						position = line.split()[1]
						samples = line.split()[9:12]
						out.write(scaffold)
						out.write("\t")
						out.write(str(position))
						out.write("\t")
						out.write("GEN:RD:AD:%:pval:ase")
						out_density.write(str(count))
						out_density.write(",")
						out_density.write(scaffold)
						out_density.write(",")
						out_density.write(str(position))
						for sample in samples:
							ase_valid = 0 
							genotype = sample.split(":")[0]
							rd_value = sample.split(":")[4]
							ad_value = sample.split(":")[5]
							if rd_value == '.' or ad_value == '.':
								continue

							RD = float(rd_value)
							AD = float(ad_value)
							ref_allele_fraction = (float(RD))/((float(RD))+(float(AD)))
							ref_allele_fraction_density = (float(RD))/((float(RD))+(float(AD)))
							if genotype == "0/1":
								p = float(binom_test_ase(RD,AD))
								if p < 0.05:
									if ref_allele_fraction >= 0.7 or ref_allele_fraction <= 0.3:
										ase_valid = 1
									else:
										ase_valid = 0
								else:
									ase_valid = 0
							else:
								p = "NA"
								ase_valid = "NA"
								ref_allele_fraction_density = "NA"
							out.write("\t")
							out.write(genotype)
							out.write(":")
							out.write(str(RD))
							out.write(":")
							out.write(str(AD))
							out.write(":")
							out.write(str(ref_allele_fraction))
							out.write(":")
							out.write(str(p))
							out.write(":")
							out.write(str(ase_valid))

							out_density.write(",")
							out_density.write(str(ref_allele_fraction_density))
						out.write("\n")
						out_density.write("\n")

if __name__ == '__main__':
	main()
