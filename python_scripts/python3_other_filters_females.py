#!/usr/bin/python3
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

#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="A vcf file")
parser.add_argument("outfile", type=str,
                    help="A filtered vcf file for triallelic SNPs and Ns and no genotype info")

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

    startTime = time.time()
    count_valid_SNPs = 0
    count_triallelic = 0
    count_N_SNPs = 0
    vcf_lines = 0
    bp = ["A","G","T","C"]

    vcf_reader = vcf.Reader(filename=args.infile)
    vcf_writer = vcf.Writer(open(args.outfile, "w"), vcf_reader)

    for record in vcf_reader:
        vcf_lines += 1
        sample_DP = 0
        genotypes = []

        # Remove INDEL info, only consider SNPs
        if record.is_snp: 
            # Check INDEL filter has worked. REF should be 1bp
            if len(record.REF) == 1:
                # Check REF is not N
                if record.REF not in bp:
                    count_N_SNPs += 1
                else:
                    # Check not triallelic site
                    if len(record.ALT) != 1:
                        count_triallelic += 1
                    else:
                        for sample in record.samples:
                            # Ignore missing data (./.)
                            if sample['GT'] != "./.":
                                # Ignore sites with quality < 20
                                if sample['GT'] == "0/1":
                                    if sample['ABQ'] >= 20 and sample['RBQ'] >= 20:
                                        genotypes.append(sample['GT'])
                                elif sample['GT'] == "0/0":
                                    if sample['RBQ'] >= 20:
                                        genotypes.append(sample['GT'])
                                elif sample['GT'] == "1/1":
                                    if sample['ABQ'] >= 20:
                                        genotypes.append(sample['GT'])
        
        if len(genotypes) == 3:
            count_valid_SNPs += 1
            vcf_writer.write_record(record)

    print("No of lines in vcf file =", vcf_lines)
    print("No of valid SNPs =", count_valid_SNPs) 
    print("No of triallelic+ SNPs =", count_triallelic)
    print("No of SNPS where REF=N =", count_N_SNPs)

    print('The script took {0} second(s)!'.format(time.time() - startTime))

if __name__ == '__main__':
    main()
