#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
'''
'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import numpy as np
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("rpkm_all_chrms", type=str,
                    help="File containing gene rpkm values for all samples")
parser.add_argument("gene_pos", type=str,
                    help="File with genes and their position in the genome")
parser.add_argument("outfile_rpkm__fc_pos", type=str,
                    help="An outfile of genes on the sex chromosome, thier position in the chromosome and their FC")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

  gene_position = defaultdict(list)

  with open(args.gene_pos, "r") as gene_pos:
    next(gene_pos)
    for line in gene_pos:
      gene = line.split(",")[1]
      position_in_genome = line.split(",")[3].rstrip()
      chromosome = line.split(",")[2]
      gene_position[gene].append(chromosome)
      gene_position[gene].append(position_in_genome)

  with open(args.outfile_rpkm__fc_pos, "w") as out:
    out.write("gene,chromosome,position_in_genome,avg_male_expr,log_male_expr,avg_female_expr,log_female_expr,fc\n")
    with open(args.rpkm_all_chrms, "r") as rpkm:
      next(rpkm)
      for line in rpkm:
        line = line.rstrip()
        gene = line.split("\t")[0]
        if gene in gene_position:
          female1 = float(line.split("\t")[1])
          female2 = float(line.split("\t")[2])
          female3 = float(line.split("\t")[3])
          male1 = float(line.split("\t")[4])
          male2 = float(line.split("\t")[5])
          male3 = float(line.split("\t")[6])

          average_male = float((male1 + male2 + male3)/3) + 1
          average_female = float((female1 + female2 + female3)/3) + 1

          log_male = float(np.log2(average_male))
          log_female = float(np.log2(average_female))

          fc = float(log_male - log_female)

          out.write(gene)
          out.write(",")
          out.write("LG"+gene_position[gene][0])
          out.write(",")
          out.write(gene_position[gene][1])
          out.write(",")
          out.write(str(average_male))
          out.write(",")
          out.write(str(log_male))
          out.write(",")
          out.write(str(average_female))
          out.write(",")
          out.write(str(log_female))
          out.write(",")
          out.write(str(fc))
          out.write("\n")



if __name__ == '__main__':
    main()
