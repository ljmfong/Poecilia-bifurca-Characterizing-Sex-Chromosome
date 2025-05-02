#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' FILTER Assembly file to remove ncRNA
Takes a folder of blast tophit files and filters the Assembly to remove these transcripts.
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
parser.add_argument("blastfolder", type=str,
                    help="Infolder containing blast output files")
parser.add_argument("assembly", type=str,
                    help="Fasta file")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".tophits")]

def get_tophits(file):
    '''Returns a list of all query transcripts in blast output file'''
    transcriptlist = []
    with open(file, "r") as infile:
        for line in infile:
            line = line.rstrip()
            transcript = line.split()[0]
            transcriptlist.append(transcript)
    return transcriptlist

def get_gene_transcript_id(file):
    gene_isoform_dict = {}
    with open(file,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith(">"):
                    gene = line.split(">")[1].rstrip()
            else:
                sequence = line.rstrip()

                gene_isoform_dict[gene] = sequence
    return gene_isoform_dict

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    #read blast files
    blast_tophits = list_folder(args.blastfolder)
    print "Number of blast tophit files =", len(blast_tophits)

    #get list of isoforms with blast hits to ncrna
    removealltranscriptlist = []
    for blast in blast_tophits:
        print blast
        removetranscriptlist = get_tophits(blast)
        removealltranscriptlist = removealltranscriptlist+removetranscriptlist
        print "No. of transcript with hits to ncrna =", len(removetranscriptlist)
    removealltranscriptlist = set(removealltranscriptlist)
    print "Total no. of transcript to remove from GTF =", len(removealltranscriptlist)
    print removealltranscriptlist

    # get transcript and sequence
    gene_sequence_dict = get_gene_transcript_id(args.assembly)
    print "Total no. of transcripts in assembly =", len(gene_sequence_dict)

    #filter and print assembly file
    count = 0
    original = 0
    outfilename = args.assembly[:-3]+"_ncrnafiltered.fa"
    print "Writing to ... ", outfilename
    with open(outfilename,"w") as outfile:
        for gene in gene_sequence_dict:
            original += 1
            if gene not in removealltranscriptlist:
                count += 1
                outfile.write(">")
                outfile.write(gene)
                outfile.write("\n")
                outfile.write(gene_sequence_dict[gene])
                outfile.write("\n")

    print "Original number of lines in assembly =", original
    print "Number of filtered lines in assembly =", count
    
if __name__ == '__main__':
    main()
