#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
"This script creates a file with the start and stop positions of genes in scaffolds"
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
from collections import OrderedDict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("gene_position", type=str,
                    help="File produced by script 10.extract-counts-annotated.py from Expression analysis")
parser.add_argument("ncrnafiltered_gtf", type=str,
                    help="File produced by scropt 06.filter-GTF-ncrna.py from Expression analysis")
parser.add_argument("coordinates", type=str,
                    help="Output file with gene coordinates in scaffolds")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if not f.endswith(".DS_Store")]


#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    genes_with_position = {}
    coord_dict = defaultdict(list)

    # Get gene-scaffold association
    with open(args.gene_position, "r") as gene_pos:
        next(gene_pos)
        for line in gene_pos:
            line = line.rstrip()
            gene = line.split(",")[0]
            scaffold = line.split(",")[1]
            genes_with_position[gene] = scaffold

    # Get start and end positions of genes in scaffolds
    with open(args.ncrnafiltered_gtf, "r") as gtf:
        for line in gtf:
            if not line.startswith("#"): 
                if line.split("\t")[2] == "transcript":
                    scaffold = line.split("\t")[0]
                    gene = line.split("\t")[8].split('gene_id "')[1].split('";')[0]
                    if gene in genes_with_position:
                        if scaffold == genes_with_position[gene]:
                            start = float(line.split("\t")[3])
                            end = float(line.split("\t")[4])
                            if start > end:
                                print "ERROR - start position bigger than end position"
                                print gene, scaffold, start, end
                                sys.exit()
                            else:
                                coord_dict[gene].append(start)
                                coord_dict[gene].append(end)
                        else:
                            print "ERROR - scaffold in gtf file and position file don't match"
                            print gene, scaffold
                            sys.exit()

    # Write output file
    with open(args.coordinates, "w") as coord:
        coord.write("scaffold\tgene,start,stop\n")
        for gene in coord_dict:
            coord.write(genes_with_position[gene])
            coord.write("\t")
            coord.write(gene)
            coord.write(",")
            coord.write(str(coord_dict[gene][0]))
            coord.write(",")
            coord.write(str(coord_dict[gene][1]))
            coord.write("\n")


if __name__ == '__main__':
    main()