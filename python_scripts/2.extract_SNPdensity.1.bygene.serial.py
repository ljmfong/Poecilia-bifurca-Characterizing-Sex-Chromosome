#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import numpy as np
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("SNP", type=str,
                    help="A file of scaffolds and SNPs")
parser.add_argument("outfilepath", type=str,
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if not f.endswith(".DS_Store")]

def list_files(current_dir):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if not name.endswith(".DS_Store"):
                f = os.path.join(path, name)
                file_list.append(f)
    return file_list

# def get_file_dict(infolders, popdict):
def get_file_dict(infolders):
    filedictionary = {}
    for infolder in infolders:
        name = os.path.basename(infolder)
        # print name
        infiles = list_files(infolder)
        for infile in infiles:
            # print infile
            if infile.endswith("stringtie.bygene"):
                filedictionary[name] = infile
    print "No. of samples =", len(filedictionary)
    # popfiledictionary = defaultdict(list)
    # for pop in popdict:
    #     for name in popdict[pop]:
    #         popfiledictionary[pop].append(filedictionary[name])
    # print "No. of populations =", len(popfiledictionary)
    # return popfiledictionary
    return filedictionary

def extract_depth(source, popsnpdict, genescaffolddict):
    count = 0
    #name, sites, snps
    with open(source, "r") as infile:
        for line in infile:
            count += 1
            line = line.rstrip()
            line = line.split("\t")
            scaffold = line[0]
            gene = line[1]
            sites = float(line[2])
            snps = float(line[3])
            row = [sites, snps]
            popsnpdict[gene].append(row)
            if gene in genescaffolddict:
                if scaffold != genescaffolddict[gene]:
                    print "error"
            else:
                genescaffolddict[gene] = scaffold
    return popsnpdict, genescaffolddict

def average_density(popsnpdict):
    #filter low coverage scaffolds
    #therefore can just loop through one dictionary
    averagedict = defaultdict(list)
    for gene in popsnpdict:
        # print gene
        if len(popsnpdict[gene]) == 3:
            count = 0
            sumsites = 0
            sumsnps = 0
            # print gene
            for sample in popsnpdict[gene]:
                # print sample
                sites = float(sample[0])
                snps = float(sample[1])
                sumsites += sites
                sumsnps += snps
                if sites > 0:
                    count += 1

            # #filter for low coverage scaffolds
            # if sumsites == 0:
            if count != 3:
                pass
            else:
                averagedensity = sumsnps/sumsites
                logaveragedensity = np.log2(averagedensity+1)
                averagedict[gene].append(sumsites)
                averagedict[gene].append(sumsnps)
                averagedict[gene].append(averagedensity)
                averagedict[gene].append(logaveragedensity)
    return averagedict
           
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    #extract snp density
    # popdict = {}

    # popdict["wingei"] = [""]
    
    # popdict["AH"] = ["J001","J003","J004","J005"]
    # popdict["AL"] = ["J032","J033","J034","J035"]

    # popdict["YL"] = ["J012","J013","J014","J015"]
    # popdict["YH"] = ["J041","J042","J043","J044"]

    # popdict["QH"] = ["J016","J017","J019","J020"]
    # popdict["QL"] = ["J051","J052","J053","J054"]

    #get dictionary of files
    SNPfolders = list_folder(args.SNP)
    print "No. of infolders =", len(SNPfolders)
    filedictionary = get_file_dict(SNPfolders)
    # print filedictionary

    #get SNP density for each sample

    snpdicty = defaultdict(list)
    genescaffolddict = {}

    for name in filedictionary:
        # print name
        # print filedictionary[name]
        snpdict, genescaffolddict = extract_depth(filedictionary[name], snpdicty, genescaffolddict)
    print "Number of genes =", len(snpdict), len(genescaffolddict)

    # calculate average SNP density
    average = average_density(snpdict)
    # print average

    #print average
    outfilename = args.outfilepath+"/Femalegenome.males.genesgtffiltered.stringtie.1.bygene"
    with open(outfilename, "w") as outfile:
        header = "Scaffold, Gene, Sumsites, Sumsnps, Snpdensity, LogSnpdensity"
        outfile.write(header)
        outfile.write("\n")
        for gene in average:
            scaffold = genescaffolddict[gene]
            outfile.write(scaffold)
            outfile.write(",")
            outfile.write(gene)
            outfile.write(",")
            outfile.write(str(average[gene][0]))
            outfile.write(",")
            outfile.write(str(average[gene][1]))
            outfile.write(",")
            outfile.write(str(average[gene][2]))
            outfile.write(",")
            outfile.write(str(average[gene][3]))
            outfile.write("\n")


    # #calculate averages
    # for pop in popfiledictionary:
    #     #get SNP density for each sample for each population
    #     popsnpdict = defaultdict(list)
    #     genescaffolddict = {}
    #     print pop, len(popfiledictionary[pop])
    #     for file in popfiledictionary[pop]:
    #         snpdict, genescaffolddict = extract_depth(file, popsnpdict, genescaffolddict)
    #     print "Number of genes =", len(snpdict), len(genescaffolddict)

    #     #calculate average SNP density
    #     average = average_density(popsnpdict)
    #     print "Number of filtered genes (>each 0) =", len(average)

    #     #print average
    #     outfilename = args.outfilepath+"/"+pop+"/Femalegenome_"+pop+".males.genesgtffiltered.stringtie.1.bygene"
    #     with open(outfilename, "w") as outfile:
    #         header = "Scaffold, Gene, Sumsites, Sumsnps, Snpdensity, LogSnpdensity"
    #         outfile.write(header)
    #         outfile.write("\n")
    #         for gene in average:
    #             scaffold = genescaffolddict[gene]
    #             outfile.write(scaffold)
    #             outfile.write(",")
    #             outfile.write(gene)
    #             outfile.write(",")
    #             outfile.write(str(average[gene][0]))
    #             outfile.write(",")
    #             outfile.write(str(average[gene][1]))
    #             outfile.write(",")
    #             outfile.write(str(average[gene][2]))
    #             outfile.write(",")
    #             outfile.write(str(average[gene][3]))
    #             outfile.write("\n")

if __name__ == '__main__':
    main()