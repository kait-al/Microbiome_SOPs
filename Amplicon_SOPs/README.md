# Burton lab pipeline for processing Illumina 16S reads
###### Kait Al, PhD. February 2023

## Overview

This SOP will take paired fastq reads and generate an ASV table with approximate taxonomy assignment. The reads must be paired, overlapping, and contain paired barcode information.
This pipeline was modified from [Greg Gloor's](https://github.com/ggloor/miseq_bin), and the [dada2 tutorial](https://benjjneb.github.io/dada2/tutorial.html).

## Get your data

London Regional Genomics Centre (LRGC) posts Illumina MiSeq reads to Basespace. Unless you have your own account, runs will be posted to the general Burton Lab account - [email me](kal@uwo.ca) for login details.

Under 'Runs' select your study, then inspect the quality of your run (Data by cycle %>=Q30). To download, go to File>Download>Run. Using the BaseSpace Sequence Hub Downloader, download the FASTQ files. They will be several GB - make sure you have enough space!

## Set up your working directory

1. Make a directory for your current study. Example: `/Users/kaital/Documents/Research/16S/Study_XYZ`
2. Make a folder called `reads` - put your .fastq.gz files from BaseSpace there.
3. Make a folder called `scripts` - put a copy of your demultiplex, dada2 workflow, and other scripts there.
4. Make a folder called `figures` - the quality profiles and error rates from dada2 will go here.
5. Fill out the `samples.txt` file according to the plate used for PCR amplification and the corresponding barcodes. An example can be found in the [example_files](/example_files). Sample IDs must be unique for each set of barcodes, and use only upper and lower case letters, numbers, or underscores. No dashes, spaces, or brackets!

## Unzip your reads file

Navigate to your reads folder in command line and run:
```
gzip -d filename_R1_001.fastq.gz
gzip -d filename_R2_001.fastq.gz
```

## Demultiplex reads

Navigate back to your main working directory. Your reads should be in the `reads` folder, and you need to make an `mmv_patterns.sh` based on your `samples.txt` file.

In R, use the following code chunk to generate fasta files with barcode info in your working directory.
These files will link the 12-mer barcodes within the reads files with their accompanying barcode number i.e. 'Golay_L1', 'Golay_R1'.
```
R
library(seqinr)

write.fasta(as.list(c('^TGCATACACTGG', '^ACTCACAGGAAT', '^GTAGGTGCTTAC', '^CAGTCGTTAAGA', '^CACTACGCTAGA', '^GCTCGAAGATTC', '^TGAACGTTGGAT', '^ATGGTTCACCCG', '^CGAGGGAAAGTC', '^TACTACGTGGCC', '^GTTCCTCCATTA', '^ACGATATGGTCA', '^TATCGACACAAG', '^AGCATGTCCCGT', '^CCAGATATAGCA', '^GTGTCCGGATTC', '^ATCGCACAGTAA', '^CAGCTCATCAGC', '^GCATATGCACTG', '^TGTAGGTGTGCT', '^ACGAGACTGATT', '^CATCAGTACGCC', '^GTATCTGCGCGT', '^TGCGTCAGCTAC')),
            c('Golay_L1', 'Golay_L2', 'Golay_L3', 'Golay_L4', 'Golay_L5', 'Golay_L6', 'Golay_L7', 'Golay_L8', 'Golay_L9', 'Golay_L10', 'Golay_L11', 'Golay_L12', 'Golay_L13', 'Golay_L14', 'Golay_L15', 'Golay_L16', 'Golay_L17', 'Golay_L18', 'Golay_L19', 'Golay_L20', 'Golay_L21', 'Golay_L22', 'Golay_L23', 'Golay_L24'),
            "fwd_barcodes.fasta", open = "w", as.string = TRUE)

write.fasta(as.list(c('^CGAGGGAAAGTC', '^TACTACGTGGCC', '^GTTCCTCCATTA', '^ACGATATGGTCA', '^TATCGACACAAG', '^AGCATGTCCCGT', '^CCAGATATAGCA', '^GTGTCCGGATTC', '^ATCGCACAGTAA', '^CAGCTCATCAGC', '^GCATATGCACTG', '^TGTAGGTGTGCT', '^ACGAGACTGATT', '^CATCAGTACGCC', '^GTATCTGCGCGT', '^TGCGTCAGCTAC', '^GTAGATCGTGTA', '^CAGCTGGTTCAA', '^AGCTGATAGTTG', '^TCTACGGCACGT', '^GCATAAACGACT', '^AAGGCGCTCCTT', '^TGCCTAAGATCG', '^CTTAGCTACTCT')),
            c('Golay_R1', 'Golay_R2', 'Golay_R3', 'Golay_R4', 'Golay_R5', 'Golay_R6', 'Golay_R7', 'Golay_R8', 'Golay_R9', 'Golay_R10', 'Golay_R11', 'Golay_R12', 'Golay_R13', 'Golay_R14', 'Golay_R15', 'Golay_R16', 'Golay_R17', 'Golay_R18', 'Golay_R19', 'Golay_R20', 'Golay_R21', 'Golay_R22', 'Golay_R23', 'Golay_R24'),
            "rev_barcodes.fasta", open = "w", as.string = TRUE)

q()
```
Proceed through the remainder of the `demultiplex_cutadapt.sh` pipeline.
Once trimmed reads have been demultiplexed and primers have been removed, navigate to the final output directory `working_directory/demultiplex_cutadapt/trimmed_reads`. From here, run `mmv_patterns.sh`.

```
#give yourself permission first if required
chmod +x mmv_patterns.sh
bash mmv_patterns.sh
```
Now your reads are demultiplexed and should be renamed to your sample IDs. You're ready to run the dada2 workflow.

## Run dada2

You should have a working copy of `dada2_workflow.R` where you've made the necessary edits based on your own data. I strongly suggest running this workflow line-by-line to ensure each step reaches a logical output.

Download the most recent version of the silva database [here](https://zenodo.org/record/4587955#.Y-7XNS370bk). The training set and species assignment should be in your working directory.

Final output will be 3 files in a new folder called `data`:
1. `cutadapt_counts.txt` this is your ASV counts table
2. `cutadapt_tax.txt` this is the taxonomy annotated to each ASV, and the sequence itself.
3. `cutadapt_track.txt` this file provides sample-wise read counts from the original input to data2 through filtering, denoising, removing chimeras, all the way to the final output.

## Cleanup

Original .fastq files from Basespace in `reads` should be compressed back to .fastq.gz:
```
gzip reads/*.fastq
```
When publishing results from sequencing data, you need to deposit your reads in a public repository. Your reads in `/Working_directory/demultiplex_cutadapt/trimmed_reads` should be uploaded to NCBI SRA prior to submitting for publication. 
All other intermediate files can be deleted including `reads/reads1_trimmed_by4.fastq.gz` and `demultiplex_cutadapt/*.fastq`.
