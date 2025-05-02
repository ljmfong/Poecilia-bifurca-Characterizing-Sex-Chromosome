import argparse
import sys
import os
from collections import defaultdict
import numpy as np
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("mf_depth_file", type=str,
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
        return "LG2"
    elif scaffold == "Scaffold_2_RagTag":
        return "LG3"
    elif scaffold == "Scaffold_3_RagTag":
        return "LG9"
    elif scaffold == "Scaffold_4_RagTag":
        return "LG5"
    elif scaffold == "Scaffold_5_RagTag":
        return "LG16"
    elif scaffold == "Scaffold_6_RagTag":
        return "LG7"
    elif scaffold == "Scaffold_7_RagTag":
        return "LG13"
    elif scaffold == "Scaffold_8_RagTag":
        return "LG1"
    elif scaffold == "Scaffold_9_RagTag":
        return "LG10"
    elif scaffold == "Scaffold_10_RagTag":
        return "LG4"
    elif scaffold == "Scaffold_11_RagTag":
        return "LG17"
    elif scaffold == "Scaffold_12_RagTag":
        return "LG6"
    elif scaffold == "Scaffold_13_RagTag":
        return "LG15"
    elif scaffold == "Scaffold_14_RagTag":
        return "LG12"
    elif scaffold == "Scaffold_15_RagTag":
        return "LG14"
    elif scaffold == "Scaffold_16_RagTag":
        return "LG11"
    elif scaffold == "Scaffold_17_RagTag":
        return "LG18"
    elif scaffold == "Scaffold_18_RagTag":
        return "LG22"
    elif scaffold == "Scaffold_19_RagTag":
        return "LG20"
    elif scaffold == "Scaffold_20_RagTag":
        return "LG8"
    elif scaffold == "Scaffold_21_RagTag":
        return "LG21"
    elif scaffold == "Scaffold_22_RagTag":
        return "LG19"
    elif scaffold == "Scaffold_23_RagTag":
        return "LG23"
    else:
        return scaffold

 
#==============================================================================
#Main==========================================================================
#==============================================================================

def main():

    coord_dict = defaultdict(list)

    with open(args.outfile, "w") as coord:
        coord.write("CHROM\tSTART\tEND\tmale_average\tfemale_average\tM_F_ratio\n")
        with open(args.mf_depth_file, "r") as mfdepth:
            next(mfdepth)
            for line in mfdepth:
                scaffold = line.split("\t")[0]
                start = line.split("\t")[1]
                end = line.split("\t")[2]
                male_average = line.split("\t")[3]
                female_average = line.split("\t")[4]
                MF_ratio = line.split("\t")[5]
                name = get_name_from_scaffold(scaffold)
    # Append the values to the dictionary
                coord_dict[MF_ratio].extend([scaffold, start, end, male_average, female_average, MF_ratio, name,])
                coord.write(name)
                coord.write("\t")
                coord.write(start)
                coord.write("\t")
                coord.write(end)
                coord.write("\t")
                coord.write(male_average)
                coord.write("\t")
                coord.write(female_average)
                coord.write("\t")
                coord.write(MF_ratio)



if __name__ == '__main__':
    main()
