library(zCompositions) # CZM

library(dplyr) # Pipe

library(ggplot2) # Plotting

library(vegan) # envfit and diversity

library(ggsci) # Colours

library(ALDEx2) # ALDEx2

#### Load data ####

counts <- read.table("data/cutadapt_counts.txt", 
                     header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, 
                     quote = "", stringsAsFactors = FALSE)
tax <- read.table("data/cutadapt_tax.txt", 
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, 
                  quote = "", stringsAsFactors = FALSE)

#### #load meta data table

metadata<-read.table("metadata.txt", header=T, row.names = 1, sep='\t', comment.char = "")

#### Filter ####

# Summarize input counts table
sum(counts) ; dim(counts) ; summary(colSums(counts)) ; summary(rowSums(counts))

# Generate a relative abundance table and remove SVs accounting for < 1% of reads in every sample
# If you care about very rare things you could lower this to 0.1% (0.001)
props <- apply(counts, 1, function(x) {x/sum(x)})
filt_by_props <- counts[, apply(props, 1, max) >= 0.01] 

# Summarize output filtered counts table
sum(filt_by_props) ; dim(filt_by_props) ; summary(colSums(filt_by_props)) ; summary(rowSums(filt_by_props))

# Create filtered data frames
filtered_counts <- counts[rownames(filt_by_props), colnames(filt_by_props)]
filtered_tax <- tax[colnames(filt_by_props), ]
filtered_meta <- metadata[rownames(filt_by_props), ]

# Can't have zeroes for downstream steps, impute them to something else logically
# Handle this problem using cmultRepl in zCompositions package
# Bayesian-Multiplicative replacement of count zeros
# method="CZM" uses multiplicative simple replacement (multRepl) on the matrix of estimated probabilities
# samples as rows 
czm <- cmultRepl(filtered_counts, label = 0, method = "CZM")

#CLR transform the data
clr <- t(apply(czm, 1, function(x) {log(x) - mean(log(x))} ))

# The output will have samples as ROWS
# Samples must be ROWs and features/OTUs as COLUMNS
# base R function to perform principal component analysis
pca <- prcomp(clr)

d.mvar <- sum(pca$sdev^2)

# Calculate the PC1 and PC2 variance
PC1 <- paste("PC1: ", round(sum(pca$sdev[1]^2)/d.mvar, 3))
PC2 <- paste("PC2: ", round(sum(pca$sdev[2]^2)/d.mvar, 3))
biplot(pca, var.axes=T, scale=0, xlab=PC1, ylab=PC2, cex=c(0.6, 0.6))
