''' PICK BEST ISOFROM
This script takes a folder of individual sample folders containing isoforms.results files from RSEM. Picks the best isoform 
for each gene. Best isoform is the isoform with the highest average expression but if there are multiple most highly 
expressed then the longest isoform is chosen. Prints list of each gene and best isoform into one file.
'''
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
parser.add_argument("infolder", type=str,
                    help="Infolder of folders containing isoforms.results files from RSEM")
parser.add_argument("reference_fasta", type=str,
                    help="A fasta file of reference sequences")
parser.add_argument("outfile", type=str,
                    help="An outfile")
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

def list_file(current_dir):
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if not name.endswith(".DS_Store"):
                if name.endswith(".isoforms.results"):
                    f = os.path.join(path, name)
                    return f

def read_scaffolds(fasta):
    '''Returns a dictionary of isoforms and sequences'''
    SG={}
    try:
        with open(fasta, "r") as infile:
            for line in infile.readlines():
                if line[0] == ">":
                    name = line[1:].rstrip().split()[0]
                    header = line[1:].rstrip()
                    SG[name] = ["", header]
                else:
                    SG[name][0] += line.rstrip()
            return SG
    except IOError:
        print "!----ERROR----!"
        print "File %s does not exit!" % fasta
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)

def get_all_length(fasta_dict):
    '''Returns dictionary of isoforms and their sequence length'''
    dict = {}
    for fasta in fasta_dict:
        dict[fasta] = float(len(fasta_dict[fasta][0]))
    return dict
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    #get isoform lengths
    fasta_dict = read_scaffolds(args.reference_fasta)
    all_length_dict = get_all_length(fasta_dict)

    #extract isoform expression
    print "Extracting expression of all isoforms ..."
    header = []
    folders = list_folder(args.infolder)
    isoform_count = defaultdict(list)
    gene_isoform_index = defaultdict(list)
    isoform_gene_index = defaultdict(list)
    for folder in folders:
        file = list_file(folder)
        if not file is None:
            print file
            sample = os.path.basename(folder)
            header.append(sample)
            with open(file, "r") as infile:
                for line in infile:
                    if not line.startswith("transcript_id"):
                        line = line.rstrip()
                        isoform = line.split("\t")[0]
                        gene = line.split("\t")[1]
                        fpkm = line.split("\t")[6]
                        isoform_count[isoform].append(fpkm)  
                        if gene in gene_isoform_index:
                            if isoform not in gene_isoform_index[gene]:
                                gene_isoform_index[gene].append(isoform) 
                        else:
                            gene_isoform_index[gene].append(isoform) 
                        if isoform in isoform_gene_index:
                            if gene not in isoform_gene_index[isoform]:
                                isoform_gene_index[isoform].append(gene) 
                        else:
                            isoform_gene_index[isoform].append(gene) 

    print "Header = ", header
    print "Total no. of isoforms =", len(isoform_count)
    print "Total no. of genes =", len(gene_isoform_index)

    #calculate averge fpkm across all samples
    print "Calculating average FPKM across samples ..."
    isoform_average = {}
    for isoform in isoform_count:
        averagelist = []
        for fpkm in isoform_count[isoform]:
            fpkm = float(fpkm)
            averagelist.append(fpkm)
        average = sum(averagelist)/float(len(averagelist))
        isoform_average[isoform] = average   

    #for each gene find the best isoform
    #best isoform is the isoform with the highest expression, or if there are multiple isoforms with the highest expr then the longest one is picked
    print "Finding best isoforms ..."
    best_isoform_dict = {}
    for gene in gene_isoform_index:
        if len(gene_isoform_index[gene]) == 1:
            #if only one isoform for a given gene, pick that isoform
            isoform = gene_isoform_index[gene][0]
            bestisoform = isoform
            best_isoform_dict[gene] = bestisoform
        else:
            #if multiple isoforms for a given gene, pick the most highly expressed
            #calculate average expression for each isoform
            average_fpkm_dict = {}
            isoforms = gene_isoform_index[gene]
            for isoform in isoforms:
                average_fpkm = float(isoform_average[isoform])
                average_fpkm_dict[isoform] = average_fpkm
            
            #order the isoforms by average expression
            average_fpkm_dict_sorted = OrderedDict(sorted(average_fpkm_dict.items(), key=lambda t: t[1], reverse=True))
            average_fpkm_sorted = list(average_fpkm_dict_sorted.values())
            highest_fpkm = float(average_fpkm_sorted[0])
            
            #count if multiple isoforms have the same highest average expression
            equal_highest_fpkm_count = sum(i >= highest_fpkm for i in average_fpkm_sorted)
            if equal_highest_fpkm_count == 1:
                #if one isoform has a higher expression that the others, pick that isoform
                bestisoform = average_fpkm_dict_sorted.items()[0][0]
                best_isoform_dict[gene] = bestisoform
            else:
                #pick the longest isoform of the most equally highly expressed isoforms
                equal_highest_fpkm_isoform_names = list(average_fpkm_dict_sorted)[:equal_highest_fpkm_count]
                length_dict = {}
                for isoform in equal_highest_fpkm_isoform_names:
                    length = float(all_length_dict[isoform])
                    length_dict[isoform] = length

                #order the most equally highly expressed isoforms by length
                length_dict_sorted = OrderedDict(sorted(length_dict.items(), key=lambda t: t[1], reverse=True))
                length_sorted = list(length_dict_sorted.values())
                longest_length = float(length_sorted[0])

                #count if multiple most highly expressed isoforms have the same length
                equal_longest_count = sum(i >= longest_length for i in length_sorted)
                if equal_longest_count == 1:
                    #pick the single longest most highly expressed isoform
                    bestisoform = length_dict_sorted.items()[0][0]
                    best_isoform_dict[gene] = bestisoform
                else:
                    print "Cannot choose isoforms for this gene"
                    print "isoforms =", gene_isoform_index[gene]
                    print "average fpkm =", average_fpkm_dict_sorted
                    print "length =", length_dict_sorted

    print "Number of genes with bestisoforms =", len(best_isoform_dict)
    with open(args.outfile, "w") as outfile:
        for gene in best_isoform_dict:
            outfile.write(gene)
            outfile.write("\t")
            outfile.write(best_isoform_dict[gene])
            outfile.write("\n")


if __name__ == "__main__":
    main()
