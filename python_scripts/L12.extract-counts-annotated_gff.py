#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' EXTRACT counts annotated
Takes a file of read counts across samples extracted from HTseq-count. Extracts read counts for genes on scaffolds
assigned to chromosomes and prints read counts into one file. Extracts positional information for genes on scaffolds
assigned to chromosomes and prints positional information into another file.
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
parser.add_argument("read_counts", type=str,
                    help="A txt file of read counts across samples")

parser.add_argument("gff", type=str,
                    help="GfF file")

parser.add_argument("gene_position_outfile", type=str,
                    help="An outfile of the position of genes within the genome")
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

    # Get coordinates for each gene from the GFF file
    gene_index = {}
    gene_position = {}

    with open(args.gff, "r") as infile:
        for line in infile:
            if not line.startswith("#"):  # Skip comment lines
                columns = line.strip().split("\t")
                feature_type = columns[2]
                attributes = columns[8]
                
                if "RNA" in feature_type:  # Focus on RNA features
                    # Parse the attributes field for the gene ID
                    attributes_dict = {
                        key_value.split("=")[0]: key_value.split("=")[1]
                        for key_value in attributes.split(";") if "=" in key_value
                    }
                    gene = attributes_dict.get("ID")  # Get the 'ID' attribute
                    
                    if gene:  # Ensure gene ID exists
                        chromosome = columns[0]
                        start = float(columns[3])
                        
                        # Update gene_index and gene_position
                        gene_index[gene] = chromosome
                        if gene in gene_position:
                            if start < gene_position[gene]:
                                gene_position[gene] = start
                        else:
                            gene_position[gene] = start

    print("Number of genes in GFF file =", len(gene_index), len(gene_position))

    #extract positional information for genes on scaffolds assigned to chromosomes
    count = 0
    #print("Keys in gene_index:", list(gene_index.keys())[:10])  # Print the first 10 keys

    with open(args.gene_position_outfile, "w") as outfile_start:
        outfile_start.write("Geneid,Chromosome,Start")
        outfile_start.write("\n")
        with open(args.read_counts, "r") as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith("ANN"):
                    gene = line.split("\t")[0]
                    print("Gene from read counts:", gene)
                    chromosome = gene_index[gene]
                    genestart = gene_position[gene]
                    #if chromosome in gene_index:
                    count += 1
                    outfile_start.write(gene)
                    outfile_start.write(",")
                    outfile_start.write(chromosome)
                    outfile_start.write(",")
                    outfile_start.write(str(genestart))
                    outfile_start.write("\n")
    print "Number of genes with count data on annotated chromosomes =", count

if __name__ == '__main__':
    main()
