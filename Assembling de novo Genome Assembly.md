## Assembling the _de novo_ Genome Assembly

This pipeline will build the _de novo_ genome assembly from the PacBio Sequencing and Illumina DNA-sequencing.

------------------------------------------------------------------------------------------------------------------------------------
###  a. PacBio:

The _de novo_ genome assembly is built using **[Canu 2.2](https://canu.readthedocs.io/en/latest/)**. You first have to convert the .bam to .fasta or .fastq file. This will output a zipped file, e.g. female_bifurca_pacbio.fastq.gz:

    venv_pacbio/bin/bam2fastq -o female_bifurca ~/bifurca_project/PacBio/demultiplex.bc2094--bc2094.bam

When running canu, note that the genome size m is representating basepairs, e.g. genomeSize=0.03m would mean 30,000 bp, and is approximate. Genome size approximately for _P. bifurca_ is 750 Mb (750,000,000bp), therefore genomSize=750.

    canu-2.2/bin/canu -p bifurca_female -d ~/bifurca_project/assembly/ genomeSize=750m corThreads=12 useGrid=false maxMemory=128 -pacbio ~/bifurca_project/PacBio/female_bifurca.fastq.gz 2>canu_2.2_assemble.txt

You can check the completeness of the initial assembly with **[BUSCO](https://busco.ezlab.org/)**. We used the _cyprinodontiformes_odb10_ lineage dataset. Something to note about the Canu assembly outputs:

    prefix.correctedReads.fasta.gz: The reads after correction.
    prefix.trimmedReads.fasta.gz: The corrected reads after overlap based trimming.
    prefix.contigs.fasta: Everything which could be assembled and is the full assembly, including both unique, repetitive, and bubble elements.
    prefix.unassembled.fasta: Reads and low-coverage contigs which could not be incorporated into the primary assembly.

To run BUSCO, example command: 

    busco -i bifurca_female.contigs.fasta -l cyprinodontiformes_odb10 -o canu_busco -m genome -f

###  b. Improving Canu Assembly with short-read DNA-sequencing

####  i. Clean-up the short-read DNA-sequencing
To improve the Canu Assembly, you can use the Illumina short-read DNAsequence. To do so, first check the quality of the intial DNA-sequencing, trim and concatenate the sequencings that were run on multiple lanes, and check the quality after trimming. First make a text file with the pathway & name to all the files you want to QC (see bifurca_rna_samples.txt for an example - ensure there is a space after each line so FASTQC knows to read each file independently):

    FILE="bifurca_dna_samples.txt"; for i in $(cat "$FILE"); do /Linux/bin/fastqc $i -o dna_fastqc_untrimmed; done

You can check multiple **[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** files with **[MultiQC](https://multiqc.info/)**:

    cd rna_fastqc_untrimmed/
    /Linux/bin/multiqc ./

Then trim the raw reads using the python script:

    python 0.trimmomatic_bifurca_DNA.py ~/bifurca_project/dna_seq/ Trimmomatic-0.36/trimmomatic-0.36.jar Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa >running_trimmomatic_DNA.txt

The output will look like this:

    Picked up _JAVA_OPTIONS: -Xms32g -Xmx32g -Djava.io.tmpdir=/home/ljmfong/tmp
    TrimmomaticPE: Started with arguments: -phred33 /home/ljmfong/bifurca_project/dna_seq/bif_DNA_female_2_DNA_S32_L003_R1_001.fastq.gz /home/ljmfong/bifurca_project/dna_seq/bif_DNA_female_2_DNA_S32_L003_R2_001.fastq.gz /home/ljmfong/bifurca_project/dna_seq/bif_DNA_female_2_DNA_S32_L003_R1_001.output_forward_paired.fq.gz /home/ljmfong/bifurca_project/dna_seq/bif_DNA_female_2_DNA_S32_L003_R1_001.output_forward_unpaired.fq.gz /home/ljmfong/bifurca_project/dna_seq/bif_DNA_female_2_DNA_S32_L003_R2_001.output_reverse_paired.fq.gz /home/ljmfong/bifurca_project/dna_seq/bif_DNA_female_2_DNA_S32_L003_R2_001.output_reverse_unpaired.fq.gz
    ILLUMINACLIP:/Linux/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
    Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
    Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    Using Long Clipping Sequence: 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
    Using Long Clipping Sequence: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    ILLUMINACLIP: Using 1 prefix pairs, 4 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences

Concatenate the lanes and then do a quality check. The goal of concatenating is ensuring that order is preserved! You want to make sure all the forward lanes are pooled together, i.e. L001_R1; L002_R1; L003_R1; L004_R1 = individual_forward; i.e. L001_R2; L002_R2; L003_R2; L004_R2 = individual_reverse. An example of the command:

    cat bif_DNA_female_1_DNA_S31_L001_R1_001.output_forward_paired.fq.gz bif_DNA_female_1_DNA_S31_L002_R1_001.output_forward_paired.fq.gz bif_DNA_female_1_DNA_S31_L003_R1_001.output_forward_paired.fq.gz bif_DNA_female_1_DNA_S31_L004_R1_001.output_forward_paired.fq.gz > ~/bifurca_project/dna_concatenated_paired_reads/bif_female_1.output_forward_paired.fq.gz

Re-run the quality checks to see how the quality has improved.

#### ii. Improving the Canu Assembly

You can first align the short-read DNA-sequences to the Canu assembly using **[BWA](https://github.com/lh3/bwa)**. First make an index of the Canu Assembly, then align using bwa and convert into a bam file with **[SAMTools](http://www.htslib.org/)**, example command below:

    bwa index ~/bifurca_project/assembly/bifurca_female.contigs.fasta
    bwa mem -t 12 ~/bifurca_project/assembly/bifurca_female.contigs.fasta ~/bifurca_project/dna_concatenated_paired_reads/bif_female_1.output_forward_paired.fq.gz ~/bifurca_project/dna_concatenated_paired_reads/bif_female_1.output_reverse_paired.fq.gz | /Linux/bin/samtools sort > Female_1_aln.bam

You will use **[Pilon](https://github.com/broadinstitute/pilon)** to correct the long-read assembly.








