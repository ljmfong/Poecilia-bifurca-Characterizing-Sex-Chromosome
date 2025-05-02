''' TRIMMOMATIC
This script takes a folder containing paired-end fastq files, finds forward and reverse pairs for each sample, 
and executes commands to quality trim them with Trimmomatic
'''
#==============================================================================
import argparse
import sys
import os
from subprocess import Popen, list2cmdline
from itertools import product
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="Infolder of fastq.gz files to trim")
parser.add_argument("pathtotrimmomatic", type=str,
                    help="This is either the trimmomatic-0.39.jar file if you have it installed in your PATH or the path to trimmomatic eg /Users/iuliadarolti/Desktop/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar")
parser.add_argument("pathtoadapters", type=str,
                    help="Path to file containing adapter sequences")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def exec_in_row(cmds):
    '''Exec commands one after the other until finished'''
    if not cmds:
        return  # empty list
    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)
    for task in cmds:
        print task
        p = Popen(task, shell=True)
        p.wait()
    if done(p):
            if success(p):
                print "done!"
            else:
                fail()

def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith("fastq.gz")]
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    infiles = list_folder(args.infolder)
    print "\nNumber of infiles (fastq.gz) =", len(infiles)

    #Creates dictionary of forward and reverse pair for each sample
    filedictionary = defaultdict(list)
    for infile in infiles:
        name = os.path.basename(infile)
        uniqid = name.split("_")[2]+"_"+name.split("_")[3]+"_"+name.split("_")[6]  #depends on naming of files; the uniqid should be something that is unique to each sample
        #In this case, we are splitting up the "bif_DNA_sex_individ_DNA_id#_lane_read_#.fastq.gz" into just the sex and inidividual number, e.g. "Female1" or "Male2"
        filedictionary[uniqid].append(infile)
    print "Name of unique IDs =", uniqid #This provides an example of what the unique names might look like
    print "Number of samples =", len(filedictionary) #check that this number is half the total number of infiles printed above as each sample has two files, a forward and a reverse

    #Check that each sample has a forward and a reverse
    for sample in filedictionary:
        check = []
        for s in filedictionary[sample]:
            uniqid = s.split("_")[-2] #depends on naming of files; the uniqid here should be something uniqe to the forward or reverse files (eg R1/R2; forward/reverse)
            check.append(uniqid)
        if len(check) == 2 and str("R1") in check and str("R2") in check: #depends on naming of files
                pass
        else:
            print "ERROR - sample", sample, "missing forward or reverse file"
            break

    #Prepare trimmomatic commands
    print "Adapter file =", args.pathtoadapters
    print "Path to trimmomatic =", args.pathtotrimmomatic
    maketrimmomaticrun = []
    for sample in filedictionary:
        for s in filedictionary[sample]:
            if s.split("_")[-2] == str("R1"): #depends on naming of files
                forward = s
            elif s.split("_")[-2] == str("R2"): #depends on naming of files
                reverse = s
        output_forward_paired = forward.split(".fastq.gz")[0] + ".output_forward_paired.fq.gz"
        output_reverse_paired = reverse.split(".fastq.gz")[0] + ".output_reverse_paired.fq.gz"
        output_forward_unpaired = forward.split(".fastq.gz")[0] + ".output_forward_unpaired.fq.gz"
        output_reverse_unpaired = reverse.split(".fastq.gz")[0] + ".output_reverse_unpaired.fq.gz"
        trimmomaticcommand = ["java -jar "+args.pathtotrimmomatic+" PE -phred33"+" "+
                                forward+" "+
                                reverse+" "+
                                output_forward_paired+" "+
                                output_forward_unpaired+" "+
                                output_reverse_paired+" "+
                                output_reverse_unpaired+" "+
                                "ILLUMINACLIP:"+args.pathtoadapters+
                                ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"]
        maketrimmomaticrun.append(trimmomaticcommand)

    #Run trimmomatic
    print "Number of commands to run =", len(maketrimmomaticrun), "\n"
    exec_in_row(maketrimmomaticrun)
    print "\n"

if __name__ == "__main__":
    main()
