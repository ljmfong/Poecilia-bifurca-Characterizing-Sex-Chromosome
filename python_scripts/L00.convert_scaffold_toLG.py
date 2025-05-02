import argparse
import sys
import os
from collections import defaultdict
import numpy as np
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("gene_pos", type=str,
                    help="The file containing gene position info")

parser.add_argument("outfile", type=str,
                    help="The outfile with the converted gene position information")

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
  
def get_name_from_scaffold(scaffold):
    '''Mapping of scaffold to name'''
    if scaffold == "Scaffold_1_RagTag":
        return "2"
    elif scaffold == "Scaffold_2_RagTag":
        return "3"
    elif scaffold == "Scaffold_3_RagTag":
        return "9"
    elif scaffold == "Scaffold_4_RagTag":
        return "5"
    elif scaffold == "Scaffold_5_RagTag":
        return "16"
    elif scaffold == "Scaffold_6_RagTag":
        return "7"
    elif scaffold == "Scaffold_7_RagTag":
        return "13"
    elif scaffold == "Scaffold_8_RagTag":
        return "1"
    elif scaffold == "Scaffold_9_RagTag":
        return "10"
    elif scaffold == "Scaffold_10_RagTag":
        return "4"
    elif scaffold == "Scaffold_11_RagTag":
        return "17"
    elif scaffold == "Scaffold_12_RagTag":
        return "6"
    elif scaffold == "Scaffold_13_RagTag":
        return "15"
    elif scaffold == "Scaffold_14_RagTag":
        return "12"
    elif scaffold == "Scaffold_15_RagTag":
        return "14"
    elif scaffold == "Scaffold_16_RagTag":
        return "11"
    elif scaffold == "Scaffold_17_RagTag":
        return "18"
    elif scaffold == "Scaffold_18_RagTag":
        return "22"
    elif scaffold == "Scaffold_19_RagTag":
        return "20"
    elif scaffold == "Scaffold_20_RagTag":
        return "8"
    elif scaffold == "Scaffold_21_RagTag":
        return "21"
    elif scaffold == "Scaffold_22_RagTag":
        return "19"
    elif scaffold == "Scaffold_23_RagTag":
        return "23"
    else:
        return scaffold

 
#==============================================================================
#Main==========================================================================
#==============================================================================

def main():

    coord_dict = defaultdict(list)

    with open(args.outfile, "w") as coord:
        coord.write("scaffold,gene,LG,start\n")
        with open(args.gene_pos, "r") as gene_pos:
            next(gene_pos)
            for line in gene_pos:
                gene = line.split(",")[0]
                scaffold = line.split(",")[1]
                pos = line.split(",")[2]
                name = get_name_from_scaffold(scaffold)
    # Append the values to the dictionary
                coord_dict[gene].extend([scaffold, pos, name,])
                coord.write(scaffold)
                coord.write(",")
                coord.write(gene)
                coord.write(",")
                coord.write(name)
                coord.write(",")
                coord.write(pos)


if __name__ == '__main__':
    main()
