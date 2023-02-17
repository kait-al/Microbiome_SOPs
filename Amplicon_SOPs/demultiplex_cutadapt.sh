#!/bin/bash

#Activate cutadapt based on how you installed it
#source /Users/kaital/cutadapt-venv/bin/activate

conda activate cutadaptenv

#Look at your fastq files! To interpret them, consult Illumina's explanation [here](https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#        Step 1:   Trim 4 base pairs from 5' before demultiplexing
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#There are 4 nucleotides before the barcodes start. Remove them.
cutadapt -u 4 -U 4 -o reads/reads1_trimmed_by4.fastq -p reads/reads2_trimmed_by4.fastq reads/Burton-Morin-319_S1_L001_R1_001.fastq reads/Burton-Morin-319_S1_L001_R2_001.fastq
