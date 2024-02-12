#!/bin/bash

#Activate cutadapt based on how you installed it
#source /Users/kaital/cutadapt-venv/bin/activate

conda activate cutadaptenv

#Look at your fastq files! To interpret them, consult Illumina's explanation (https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html)
#The second line is the actual base call sequence.
#Illumina adapters are normally removed already.
#Do you have 4 nucleotides before the barcode? If so, perform step 1. If not, skip to step 2.

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#        Step 1:   Trim 4 base pairs from 5' before demultiplexing
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#There are 4 nucleotides before the barcodes start. Remove them.
cutadapt -u 4 -U 4 -o reads/reads1_trimmed_by4.fastq -p reads/reads2_trimmed_by4.fastq reads/your_filename_R1_001.fastq reads/your_filename_R2_001.fastq

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#        Step 2:   Demultiplex - for combinatorial barcodes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Adjust error(-e) as you see fit.  Cutadapt picks the best fitting sample match
# Barcodes should be anchored with ^ (i.e. ^ACGTGCCT...etc)
# Removes barcodes, renames files based on fwd and rev barcodes

mkdir demultiplex_cutadapt
cutadapt \
    -e 0.15 --no-indels --discard-untrimmed \
    -g file:fwd_barcodes.fasta \
    -G file:rev_barcodes.fasta \
    -o demultiplex_cutadapt/{name1}-{name2}.1.fastq \
    -p demultiplex_cutadapt/{name1}-{name2}.2.fastq \
    reads/reads1_trimmed_by4.fastq reads/reads2_trimmed_by4.fastq



cd demultiplex_cutadapt

mkdir trimmed_reads

# remove V4 primers
for i in *.1.fastq
do
SAMPLE=$(echo ${i} | sed "s/\.1\.fastq//")
cutadapt \
-e 0.1 -m 1 --discard-untrimmed -j 4 \
-a GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
-A GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
-o trimmed_reads/${SAMPLE}.1.fastq.gz \
-p trimmed_reads/${SAMPLE}.2.fastq.gz \
${SAMPLE}.1.fastq ${SAMPLE}.2.fastq
done

conda deactivate
# for V3-V4 BAKT primers 
#-a CCTACGGGNGGCWGCAG...ATTAGATACCCBDGTAGTC \
#-A GACTACHVGGGTATCTAAT...CTGCWGCCNCCCGTAGG \
