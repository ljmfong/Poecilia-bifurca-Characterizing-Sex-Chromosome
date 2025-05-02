#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' EXTRACTS bestisoform sequences
Takes a reference fasta file and list of genes and bestisoforms. 
Picks sequence for the best isoform for each gene into a new file.
'''
#==============================================================================
import argparse
import sys
import csv
import os
from itertools import combinations
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("fasta", type=str,
					help="A fasta file")
parser.add_argument("bestisoforms", type=str,
					help="A text file of gene and bestisoform")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def read_scaffolds(fasta):
	SG={}
	try:
		with open(fasta, "r") as infile:
			for line in infile.readlines():
				if line[0] == ">":
					name = line[1:].rstrip().split()[0]
					header = line[1:].rstrip()
					SG[name] = ["", header]
				else:
					SG[name][0] += line.rstrip()
			return SG
	except IOError:
		print "!----ERROR----!"
		print "File %s does not exit!" % fasta
		sys.exit(1)
	except KeyboardInterrupt:
		sys.exit(1)

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	#read in fasta file of sequences
	fasta_dict = read_scaffolds(args.fasta)
	print "Number of sequences =", len(fasta_dict)

	#read in text file of bestisoforms
	bestisoforms_dict = {}
	with open(args.bestisoforms, "r") as infile:
		for line in infile:
			line = line.rstrip()
			gene = line.split("\t")[0]
			bestisoform = line.split("\t")[1]
			bestisoforms_dict[bestisoform] = gene
	print "Number of best isoforms =", len(bestisoforms_dict)

	#print only best isoform sequences
	if args.fasta.endswith(".fa"):
		outfile = args.fasta.split(".fa")[0]+".bestisoform.fa"
	elif args.fasta.endswith(".fasta"):
		outfile = args.fasta.split(".fasta")[0]+".bestisoform.fa"
	else:
		print "Fasta file not named correctly (ie .fa or .fasta)"
		sys.exit()

	count = 0
	with open(outfile, "w") as outfile:
		for fasta in fasta_dict:
			if fasta in bestisoforms_dict:
				count += 1
				outfile.write(">"+fasta)
				outfile.write("\n")
				outfile.write(fasta_dict[fasta][0])
				outfile.write("\n")
	print "Number of best isoforms printed =", count


if __name__ == '__main__':
	main()
