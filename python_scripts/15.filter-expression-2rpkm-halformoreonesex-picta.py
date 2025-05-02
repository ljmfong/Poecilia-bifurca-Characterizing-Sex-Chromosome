#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' FILTER expression
Takes a file of genes and their rpkm values for each sample as output from edgeR.
Filters expression. Gene must be expressed > 2FPKM in half or more of the individuals.
The script creates a list of genes that have passed the filtering threshold and outputs 
a file containing genes that have passed the filtering threshold and their rpkm values. 

Takes a file of read counts across samples extracted from HTseq-count. (Header starts with "Geneid", and 
each gene name starts with "MSTRG"). Writes file with read counts for all genes in the list that have passed 
the 2rpkm filtering threshold.
'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("rpkm_file", type=str,
                    help="File containing gene rpkm values for all samples")
parser.add_argument("read_counts_file", type=str,
                    help="File containing gene read counts for all samples")
parser.add_argument("outfile_rpkm", type=str,
                    help="An outfile of filtered genes and their rpkm expression")
parser.add_argument("outfile_counts", type=str,
                    help="An outfile of filtered genes and their read counts")
parser.add_argument("conditions", type=str,
                    help="List of sex of each sample eg 'F,F,F,M,M,M' ")
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

  numberofsamples_Male = 0
  numberofsamples_Female = 0

  # count = 0
  
    
  count_all_genes = 0
  count_filtered_genes = 0
  count_filtered_check = 0

  genespassingthreshold = []
   
  with open(args.outfile_rpkm, "w") as outfile_rpkm:
      with open(args.rpkm_file, 'r') as infile:
          for line in infile:
              line = line.rstrip()
              # Printing header
              if not line.startswith("ANN"):
                  outfile_rpkm.write(line)
                  outfile_rpkm.write("\n")
                  conditions = args.conditions.split(",")
                  for sex in conditions:
                    if sex == "M":
                      numberofsamples_Male += 1
                    if sex == "F":
                      numberofsamples_Female += 1
                    # numberofsamples = (len(line.split("\t")) - 1) / 2 
                  print "number of Male samples = ", numberofsamples_Male
                  print "number of Female samples =", numberofsamples_Female
              else:
                  conditions = args.conditions.split(",")
                  
                  if line.startswith("ANN"):
                      count_all_genes += 1
                      expression_Females = []
                      expression_Males = []
                      # filtering lowly expressed genes and make list of genes passing threshold
                      expression = line.split("\t")[1:]
                      a = 0
                      total_nosamples = numberofsamples_Male + numberofsamples_Female
                      while a < total_nosamples:
                        if conditions[a] == "F":
                          expression_Females.append(expression[a])
                        else:
                          if conditions[a] == "M":
                            expression_Males.append(expression[a])
                        a += 1
                      expression_Females = [float(i) for i in expression_Females]
                      expression_Males = [float(i) for i in expression_Males]
                      numberpassingthreshold_Females = float(len([i for i in expression_Females if i > 2]))
                      numberpassingthreshold_Males = float(len([i for i in expression_Males if i > 2]))
                      if numberpassingthreshold_Females/numberofsamples_Female >= 0.5 or numberpassingthreshold_Males/numberofsamples_Male >= 0.5:
                        gene = line.split("\t")[0]
                        genespassingthreshold.append(gene)
                        count_filtered_genes += 1
                        #Write file with genes passing threshold and their rpkm values
                        outfile_rpkm.write(line)
                        outfile_rpkm.write("\n")
                          
      # create new read counts file containing genes passing threshold only
      with open(args.outfile_counts, 'w') as outfile_counts:
        with open(args.read_counts_file, 'r') as infile:
          for line in infile:
            if not line.startswith("ANN"):
              outfile_counts.write(line)
            else:
              gene = line.split("\t")[0]
              if gene in genespassingthreshold:
                count_filtered_check += 1
                outfile_counts.write(line)
          if count_filtered_check == count_filtered_genes:
            print "YES"

      print "Total number of genes = ", count_all_genes
      print "Number of genes after filtering = ", count_filtered_genes

if __name__ == '__main__':
    main()
