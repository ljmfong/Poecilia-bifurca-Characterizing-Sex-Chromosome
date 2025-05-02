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
parser.add_argument("ase_info", type=str,
					help="ase info file")
parser.add_argument("pvaladj", type=str,
					help="pvaladj file")
parser.add_argument("ase_info_adj", type=str,
					help="ase info file with adjusted pvalues")
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

	dicty = defaultdict(list)
	with open(args.pvaladj, "r") as pval:
		next(pval)
		for line in pval:
			line = line.rstrip()
			ids = line.split("\t")[1].split('"')[1]
			sample = ids.split("ase_")[1].split(":")[0]+":"+ ids.split("ase_")[1].split(":")[1]+":"+ ids.split("ase_")[1].split(":")[2]+":"+ ids.split("ase_")[1].split(":")[3]
			adj_pval = line.split("\t")[3]
			dicty[ids].append(sample)
			dicty[ids].append(adj_pval)

	count = 0
	with open(args.ase_info_adj, "w") as out:
		with open(args.ase_info, "r") as infile:
			for line in infile:
				count += 1
				if count == 1:
					out.write(line)
				else:
					line = line.rstrip()
					scaffold = line.split("\t")[0]
					position = line.split("\t")[1]
					format = line.split("\t")[2]
					s1 = line.split("\t")[3]
					s2 = line.split("\t")[4]
					s3 = line.split("\t")[5]
					out.write(scaffold)
					out.write("\t")
					out.write(str(position))
					out.write("\t")
					out.write(format)
					out.write("\t")
					for ids in dicty:
						scaffold_ids = ids.split("_")[0]
						position_ids = ids.split("_")[1]
						if scaffold == scaffold_ids and position == position_ids:
							sample_ids = ids.split("_")[-1]
							if s1 == sample_ids:
								pvaladj = float(dicty[ids][1])
								if pvaladj <= 0.05:
									switch = 1
								else:
									switch = 0
								s1 = dicty[ids][0]+":"+dicty[ids][1]+":"+str(switch)
							if s2 == sample_ids:
								pvaladj = float(dicty[ids][1])
								if pvaladj <= 0.05:
									switch = 1
								else:
									switch = 0
								s2 = dicty[ids][0]+":"+dicty[ids][1]+":"+str(switch)
							if s3 == sample_ids:
								pvaladj = float(dicty[ids][1])
								if pvaladj <= 0.05:
									switch = 1
								else:
									switch = 0
								s3 = dicty[ids][0]+":"+dicty[ids][1]+":"+str(switch)
					out.write(s1)
					out.write("\t")
					out.write(s2)
					out.write("\t")
					out.write(s3)
					out.write("\n")



 

if __name__ == '__main__':
	main()
