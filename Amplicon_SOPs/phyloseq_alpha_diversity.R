#alpha diversity metrics in phyloseq
#perform on dataset that isn't heavily filtered for max accuracy of all the metrics

library(phyloseq)

library(dplyr)

library(microbiome)

#build phyloseq object from untrimmed counts file
dm<-read.table("data/cutadapt_counts.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")

tax<-read.table("data/cutadapt_tax.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
#from dyplyr, remove sequence column
tax<-select(tax, -Sequence)
tax<-as.matrix(tax) 

OTU = otu_table(dm, taxa_are_rows = FALSE)
TAX = tax_table(tax)

meta<-read.table("metadata.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
sampledata = sample_data(meta)

physeq = phyloseq(OTU, TAX, sampledata)
physeq
# phyloseq-class experiment-level object
#
#
#

#calculate alpha diversity for each individual sample
#measures = NULL means all measures will be calculated (shannon's, chao1, simpson, etc)
div<- estimate_richness(physeq, split = TRUE, measures = NULL)

#from microbiome R package, calculate the Berger-Parker dominance index
dom<- dominance(physeq,  relative = TRUE, aggregate = FALSE)

#merge the data
div_all <- data.frame(div, dom)

# Write out the file
#import the "type" column to sort samples in graphpad
write.table(div_all, file="data/alpha_diversity_phyloseq.txt", sep="\t")

# you can go through the phyloseq tutorials with your data in physeq format
# https://joey711.github.io/phyloseq/import-data.html
