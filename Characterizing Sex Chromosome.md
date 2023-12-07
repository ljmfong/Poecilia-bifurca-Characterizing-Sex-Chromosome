## Characterizing the Sex Chromosome in _P. bifurca_

There are many different analyses that can be used to identify sex chromosomes. Here, we use a combination of male to female coverage analysis, single-nucleotide polymorphisms (SNP) density, and allele-specific expression to test for dosage compensation.

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







