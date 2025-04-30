## Building the Transcriptome and finding Allele-Specific Expression (ASE)


#### a. Checking RNA quality 
**[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**: 
Make a text file with the pathway & name to all the files you want to QC (see bifurca_rna_samples.txt for an example - ensure there is a space after each line so FASTQC knows to read each file independently): Make an output folder for your fastqcs:

    FILE="bifurca_rna_samples.txt"; for i in $(cat "$FILE"); do /Linux/bin/fastqc $i -o rna_fastqc_untrimmed; done

Example if you want to only do one file:

    /Linux/bin/fastqc rna_seq/bif_Female_1_RNA_S23_L002_R1_001.fastq.gz -o rna_fastqc_untrimmed/

Check multiple fastqc using multiqc:
**[MultiQC](https://multiqc.info/)**: Merge multiple FastQC files. Ensure that all the fastqcs are in the same directory.


    cd rna_fastqc_untrimmed/
    /Linux/bin/multiqc ./

* Important QC reports include: Per base sequence quality, Per sequence quality scores, Per seq GC content, Overrepresented sequences

#### b. RNA Trimming

**[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**: Trim off adaptors. Example command:

Using the script to trim multiple files at once: 

    /Linux/bin/python ~/bin/01.trimmomatic_bifurca_RNA.py ~/bifurca_project/rna_seq/ /Linux/Trimmomatic-0.36/trimmomatic-0.36.jar /home/ljmfong/bifurca_project/universal_adapter.fa >trimmomatic_universal.txt

Running just trimmomatic for one file at a time:

    java -jar /Linux/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 sample_1_R1.fastq.gz sample_1_R2.fastq.gz sample_1_output_forward_paired.fq.gz sample_1_output_forward_unpaired.fq.gz sample_1_output_reverse_paired.fq.gz sample_1_output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

* Run FastQC on the trimmed sequences and compare results.
    
