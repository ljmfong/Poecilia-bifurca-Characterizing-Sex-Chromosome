## Building the Transcriptome and finding Allele-Specific Expression (ASE)


#### a. Checking RNA quality and Trimming
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

