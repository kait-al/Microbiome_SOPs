library(zCompositions) # CZM
library(dplyr)
library(tibble)
library(ggplot2)

#### Load data ####

counts <- read.table("data/cutadapt_counts.txt", 
                     header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, 
                     quote = "", stringsAsFactors = FALSE)
tax <- read.table("data/cutadapt_tax.txt", 
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, 
                  quote = "", stringsAsFactors = FALSE)

#### #load meta data table

metadata<-read.table("metadata.txt", header=T, row.names = 1, sep='\t', comment.char = "")
metadata<-tibble::rownames_to_column(metadata, "SampleID")

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
d.pcx <- prcomp(clr)

d.mvar <- sum(d.pcx$sdev^2)

# Calculate the PC1 and PC2 variance
PC1 <- paste("PC1: ", round(sum(d.pcx$sdev[1]^2)/d.mvar, 3))
PC2 <- paste("PC2: ", round(sum(d.pcx$sdev[2]^2)/d.mvar, 3))

loadings<- data.frame(Variables=rownames(d.pcx$rotation), d.pcx$rotation)
values<-merge(d.pcx$x[,c(1,2)], metadata[,c("SampleID","Cohort")],
                      by.x="row.names", by.y="SampleID", all=F)
 
values$Cohort<-factor(values$Cohort, levels=c("StoneFormer","HealthyControl"))
  
theme_new <- theme_set(theme_bw())
cols <- c("StoneFormer" = "#03314B", "HealthyControl" = "#F38992")
 
ggplot(values, aes(x = PC1, y = PC2)) +
geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*80), yend = (PC2*80)),
               arrow = arrow(length = unit(5/20, "picas")),
               color = "darkgrey",
               inherit.aes = FALSE, size=0.3) +
geom_point(data = values, aes(color=Cohort)) +
#geom_text(data = values, aes(label=Row.names, hjust=-0.2), size=2)+
scale_color_manual(values = cols) +
guides(fill = guide_legend(override.aes=list(shape=21)))+
annotate("text", x = (loadings$PC1*80), y = (loadings$PC2*80), label = loadings$Variables, size=2) +
stat_ellipse(aes(x = PC1, y = PC2, colour=Cohort), data = values, geom = "path", position = "identity", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE) +
xlab(paste0("PC1: ", round(100*(d.pcx$sdev[1]^2/sum(d.pcx$sdev^2)),1),"%")) +
ylab(paste0("PC2: ", round(100*(d.pcx$sdev[2]^2/sum(d.pcx$sdev^2)),1),"%")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
