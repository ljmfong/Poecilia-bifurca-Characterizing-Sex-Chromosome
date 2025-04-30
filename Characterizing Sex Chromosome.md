## Characterizing the Sex Chromosome in _P. bifurca_

There are many different analyses that can be used to identify sex chromosomes. Here, we use a combination of male to female coverage analysis and k-mer analysis. Scripts used in this pipeline can be found under the following folder:

------------------------------------------------------------------------------------------------------------------------------------

###  a. Coverage Analysis

Following the _de novo_ genome assembly, we first aligned all the short-read DNA-sequence data to the _de novo_ genome assembly using **[BWA](https://bio-bwa.sourceforge.net/bwa.shtml)**. The following is an example command, remember to build a BWA index. Remember for paired-end data to do it for both forward and reverse sequences:

    bwa index ~/bifurca_project/ragtag_output/ragtag_scaffolds_only.fasta
    bwa aln ~/bifurca_project/ragtag_output/ragtag_scaffolds_only.fasta ~/bifurca_project/dna_concatenated_paired_reads/bif_female_1.output_forward_paired.fq.gz > Female_1_forward.sai
    bwa aln ~/bifurca_project/ragtag_output/ragtag_scaffolds_only.fasta ~/bifurca_project/dna_concatenated_paired_reads/bif_female_1.output_reverse_paired.fq.gz > Female_1_reverse.sai

Create a SAM file using bwa and find _unique_ mapping reads, the following is an example command:

    bwa sampe ~/bifurca_project/ragtag_output/ragtag_scaffolds_only.fasta ~/bifurca_project/bwa_cov_analysis/Female_1_forward.sai ~/bifurca_project/bwa_cov_analysis/Female_1_reverse.sai ~/bifurca_project/dna_concatenated_paired_reads/bif_female_1.output_forward_paired.fq.gz ~/bifurca_project/dna_concatenated_paired_reads/bif_female_1.output_reverse_paired.fq.gz > Female_1_sampe.sam
    grep 'XT:A:U' Female_1_sampe.sam > Female_1_uniq.sam

You can then estimate coverage with SOAP.coverage (see attached script). The following is an example command:

    bin/soap.coverage -sam -cvg -refsingle ~/bifurca_project/ragtag_output/ragtag_scaffolds_only.fasta -i ~/bifurca_project/bwa_sampe/Female_1_uniq.sam -o ~/bifurca_project/soapcoverage_uniq/Female_1_soapCov.txt -p 12 -window Female_1_windows.txt 50000

You can extract the coverage using the 06.extract_coverage.py script. NOTE: You need to have the soapCov.txt files in two separate directories - one for females and one for males. Also NOTE: The CovList shows the coverage of the inputted individuals (i.e. there were will be 3 values separated by "|" because I had 3 females and 3 males. In the python script, for line 28 - make sure the files in the folder match what the file name ends with - e.g. "windows.txt" or "extractcov.txt".


    python 06.extract_coverage.py ~/bifurca_project/Female_soapcov/ ~/bifurca_project/coverage_analysis/female_coverage_only.txt
    python 06.extract_coverage.py ~/bifurca_project/Male_soapcov/ ~/bifurca_project/coverage_analysis/male_coverage_only.txt

You can then manually calculate the fold change in Excel or R, and plot using the according Rscript.


###  b. Kmer Analysis

Find out how many basepairs your females and male reads contained:

First concatenate your reads if you haven't done so, example command:

    cat bif_female_1.output_forward_paired.fq.gz bif_female_1.output_reverse_paired.fq.gz >bif_female_1_both.fq.gz
    cat bif_female_1_both.fq.gz bif_female_2_both.fq.gz bif_female_3_both.fq.gz > all_bif_female_DNA.fq.gz


Then run **[Jellyfish](https://github.com/gmarcais/Jellyfish)** to get a histogram to know the profile (your expectations):

    /Linux/bin/jellyfish-2.3.3 count -C -m 31 -s 100M -t 12 <(zcat all_bif_female_DNA.fq.gz)
    /Linux/bin/jellyfish-2.3.3 count -C -m 31 -s 100M -t 12 <(zcat all_bif_male_DNA.fq.gz)

    /Linux/bin/jellyfish-2.3.3 histo ~/bifurca_project/female_mer_counts.jf > all_female_bif_counts.histo
    /Linux/bin/jellyfish-2.3.3 histo ~/bifurca_project/male_mer_counts.jf > all_male_bif_counts.histo


You can further do a k-mer comparison between species following **[HAWK](https://github.com/atifrahman/HAWK)** 

    bash countKmers 2>countKmers.log
    bash runHawk 2>runHawk_firstrun.log
    /Linux/HAWK-1.7.0-beta/supplements/EIG6.0.1-Hawk/bin/smartpca -p /Linux/HAWK-1.7.0-beta/parfile.txt
    bash post_runHawk 2>post_runHawk.log

* Edit the scripts accordingly (e.g. location and name of your files) since they are BASH scripts
* Organize the sorted_files.txt, total_kmer_counts.txt, and gwas_info.txt to represent each family and save them into each family's individual directories (runHawk is in the same directory that has the "sorted_files.txt" and "total_kmer_counts.txt"):
You will also have to make the manual case and control total_kmer_count files


At times, HAWK may not always properly produce the sums. To ensure that they have done so correctly, use the awk command to make sure they have summed correctly - example command:

    
    awk '{print $1,$5,$6,$7}' case_out_w_bonf_sorted.kmerDiff | awk '{print $1; for(i=1; i<=NF;i++) j+=$i; print j; j=0 }'  | sed 'N;s/\n/\t/' > sum_case_out_w_bonf_sorted.kmerDiff
    awk '{print $1,$8,$9,$10}' case_out_w_bonf_sorted.kmerDiff | awk '{print $1; for(i=1; i<=NF;i++) j+=$i; print j; j=0 }'  | sed 'N;s/\n/\t/' > sum_case_control_out_w_bonf_sorted.kmerDiff
    cut -f2 sum_case_control_out_w_bonf_sorted.kmerDiff | sum_case_out_w_bonf_sorted.kmerDiff - > corrected_sums_case_out_w_bonf_sorted.txt

* Note that CASE = Males, and CONTROL = Females. The header of .kmerDiff files looks like the following:

    Kmer-seq	sum_case_kmers	sum_control_kmers	p-value	case_ind_1	case_ind_2	case_ind_3	control_ind_1	control_ind_2	control_ind_3

Now you have to normalize and compare across species, since there are coverage differences between them - Example commands for 20X normalization:

    awk -F " " '{print $1,(($2/107483903941)*25000000000),(($3/100794792290)*25000000000)}' corrected_sums_case_out_w_bonf_sorted.txt >  normalized_corrected_case_w_bonf.txt
    awk '{if($3 == "0"){print;}}' normalized_corrected_case_w_bonf.txt > case_normalized_uniq.kmerDiff
    awk -F " " '($2 > 20){print;}'case_normalized_uniq.kmerDiff> case_normalized_uniq20.kmerDiff

You can check the number of unique k-mers based on the normalization 

    wc -l case_normalized_uniq20.kmerDiff

Once you have done this for all species, you can then compare the shared unique k-mers between species and the sexes - example commands:
    
    cut -f1 -d " " case_normalized_uniq20.kmerDiff > bifurca_case_normalized_uniq20_cut.kmerDiff
    sort bifurca_case_normalized_uniq20_cut.kmerDiff > bifurca_case_normalized_uniq20_cut_sorted.kmerDiff
    comm -12 picta_case_normalized_uniq20_cut_sorted.kmerDiff bifurca_case_normalized_uniq20_cut_sorted.kmerDiff > shared_PicBif_Male.txt



