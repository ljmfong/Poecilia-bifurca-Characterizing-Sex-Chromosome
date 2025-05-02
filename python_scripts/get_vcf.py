#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("consistentase", type=str,
                    help="")
parser.add_argument("vcf", type=str,
                    help="")
parser.add_argument("out", type=str,
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
def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".txt")]

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    dicty_scaff_pos = defaultdict(list)

    with open(args.consistentase, "r") as ase:
        next(ase)
        for line in ase:
            scaffold = line.split("\t")[0]
            position = float(line.split("\t")[1])
            dicty_scaff_pos[scaffold].append(position)
    # print dicty_scaff_pos
    with open(args.out, "w") as out:
        with open(args.vcf, "r") as vcf:
            for line in vcf:
                if line.startswith("#"):
                    out.write(line)
                else:
                    scaffold = line.split()[0]
                    pos = float(line.split()[1])
                    if pos in dicty_scaff_pos[scaffold]:
                        out.write(line)


if __name__ == '__main__':
    main()