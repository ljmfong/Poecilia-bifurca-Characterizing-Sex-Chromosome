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
parser.add_argument("aseinfo", type=str,
					help="")
parser.add_argument("outfile", type=str,
					help="")
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

	count_consistent_valid_sites = 0 
	with open(args.outfile, "w") as out:
		out.write("scaffold\tposition\tformat\tSample_1,2\tSample_2,2\tSample_3,2\n")
		with open(args.aseinfo, "r") as ase:
			next(ase)
			for line in ase:
				count_hetero = 0
				count_valid_ase = 0
				line = line.rstrip()
				scaffold = line.split("\t")[0]
				position = line.split("\t")[1]
				samples = line.split("\t")[3:5]
				for sample in samples:
					genotype = sample.split(":")[0]
					if genotype == "0/1":
						count_hetero += 1
						ase_valid = sample.split(":")[-1]
						if ase_valid == "1":
							count_valid_ase += 1
				if count_valid_ase !=0 and count_hetero !=0:
					if count_valid_ase == count_hetero:
						count_consistent_valid_sites += 1
						out.write(line)
						out.write("\n")
	print("Sites that show consistent ASE pattern across all samples = "), print(count_consistent_valid_sites)








if __name__ == '__main__':
	main()
