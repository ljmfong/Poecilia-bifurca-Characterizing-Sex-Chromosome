#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
#       exclude-snp-clusters.py
#==============================================================================
import argparse
import sys
import vcf
import os
# VCF requires Rpy2 internally and needs the correct R_HOME
os.environ["R_HOME"] = "/Linux/anaconda2/bin/R"
from fnmatch import fnmatch
# Make utilities folder available
sys.path.append(os.path.abspath("../"))

from multiprocessing import Process
#from utils import list_files
from collections import defaultdict

#==============================================================================
# Command line options=========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("input", type=str,
                    help="A vcf file (version 4.0 or 4.1) or folder of vcf files")
parser.add_argument("-l", "--read_length", type=int,
                    help="Length of the RNA-Seq reads")
parser.add_argument("-m", "--mismatches", type=int,
                    help="Allowed mismatches in the bowtie2 alignment")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

#==============================================================================



def find_snp_cluster(snp_positions, read_length, mismatches):
    '''Excludes clusters of neighbouring SNPs as described in Stevenson et al.,
       BMC Genomics 2013. If more SNPs occur in a window that equals the length
       of one read than the allowed number of mismatches in the alignment the
       SNPs are excluded. This should reduce the systematic bias in measures of
       ASE using RNA-Seq date with one reference genome.
    '''
    exclude = []

    for i in range(len(snp_positions)):

        cluster = [snp_positions[i]]

        for x in range(i+1, len(snp_positions)):

            if snp_positions[x] > snp_positions[i] + read_length:

                if len(cluster) > mismatches:
                    exclude += cluster
                    # print cluster
                break

            else:
                cluster.append(snp_positions[x])

    return set(exclude)


def get_snp_positions(vcf_reader):
    '''Extract the positions of SNPs from a VCF file for each chromosome.
       Returns a dictionary in the form CHR -> [p1, p2, p3..pN]'''
    snp_positions = defaultdict(list)

    for record in vcf_reader:
        snp_positions[record.CHROM].append(record.POS)
    return snp_positions


def get_snp_cluster(snp_positions, read_length, mismatches):
    '''Takes a dictionary with keys=chromosomes and list of snp positions as
       input and uses find_snp_cluster() to determine which SNPs occur in
       clusters.
    '''
    snp_cluster = defaultdict(list)

    for chromosome in snp_positions:
        snps = snp_positions[chromosome]
        snp_cluster[chromosome] = find_snp_cluster(snps,
                                                   read_length,
                                                   mismatches)

    return snp_cluster


def run_filter_clusters(args):
    def filter_clusters(infile, outfile, read_length, mismatches):
        vcf_reader = vcf.Reader(open(infile, "r"))
        snp_positions = get_snp_positions(vcf_reader)
        snp_cluster = get_snp_cluster(snp_positions, read_length, mismatches)
        vcf_reader = vcf.Reader(open(infile, "r"))
        vcf_writer = vcf.Writer(open(outfile, 'w'), vcf_reader)
        for record in vcf_reader:
            if record.POS not in snp_cluster[record.CHROM]:
                vcf_writer.write_record(record)
    processes = []
    for arg in args:
        p = Process(target=filter_clusters, args=arg)
        p.start()
        processes.append(p)
    for p in processes:
        p.join()


def handle_input(input):
    """If a directory is passed return all the files in the directory and sub-
     directories in a list, if a file is provided return a list with one entry.
    """
    if os.path.isdir(input):
        return list_files(input)
    elif os.path.isfile(input):
        return [input]


def main():

    arg_list = []

    infile = args.input

    head, tail = os.path.split(infile)
    outfile = tail.split(".")[0]+"_noclusters"+".vcf"
    outfile = os.path.join(head, outfile)
    print("Filtering file: %s" % infile)
    print(outfile)
    arg_list.append((infile, outfile, args.read_length, args.mismatches))

    # Run the filtering step in parallel
    run_filter_clusters(arg_list)


if __name__ == "__main__":
    main()
