library(ggplot2)
library(dplyr)


#### Load data ####

counts <- read.table("cutadapt_counts.txt", 
                     header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, 
                     quote = "", stringsAsFactors = FALSE)
tax <- read.table("cutadapt_tax.txt", 
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, 
                  quote = "", stringsAsFactors = FALSE)
tax <- select(tax, -Sequence)

# Generate a relative abundance table and remove SVs accounting for < 1% of reads in every sample
# If you care about very rare things you could lower this to 0.1% (0.001)
props <- apply(counts, 1, function(x) {x/sum(x)})
filt_by_props <- counts[, apply(props, 1, max) >= 0.01]


# Create filtered data frames
filtered_counts <- counts[rownames(filt_by_props), colnames(filt_by_props)]
filtered_tax <- tax[colnames(filt_by_props), ]

#transpose so samples as columns, SVs as rows
t_counts<-t(filtered_counts)

#to get rownames as genus (all same genus aggregated together)
# or aggregate at any taxonomic level in tax
rownames(t_counts)<-filtered_tax$Genus
t_counts<-as.matrix(t_counts)
t_counts_agg<-aggregate(t_counts, list(row.names(t_counts)), sum)

dim(t_counts)
# 178 108
dim(t_counts_agg)
# 60 109
#Genus name now stored in column "Group.1"

rownames(t_counts_agg) <- t_counts_agg$Group.1
t_counts_agg$Group.1 <- NULL

#If you want a specific order of samples or spaces between blocks, export this file then edit/rearrange columns in excel
#write.table(t_counts_agg, file="counts_agg_genus.txt", sep='\t', quote=F)
#add columns called "blank1", "blank2", etc. that have only 0s where you want blank blocks
#then reload the table back in as t_counts_agg
#t_counts_agg <- read.table("counts_agg_genus.txt", header=T, sep="\t", row.names=1, comment.char="", skip=0, check.names=FALSE)

#############

y <- apply(t_counts_agg, 2, function(x) { x / sum(x) } )
#1% abundance
abund <- 0.01
y2 <- y[order(rowSums(y), decreasing = TRUE),]

#check sample columns sums to 1 (100%)
colSums(y2)
dim(y2)
bugnames<-rownames(y2)

# WITHIN each sample, sum any <1% taxa into the remainder. This is for visual simplification
# Compare the plot before and after this step
y3 <- apply( y2 , 2, function(x) {
  small <- ( x < abund ) #& ( bugnames != "rem" )
  rare <- sum( x[small] )
  x[small] <- 0 # *** NA
  x["remainder"] <- x["remainder"] + rare
  return(x)
  #           return(rare) #to this to get ONLY the rare organisms
} )

pal <- colorRampPalette(colors = c("steelblue3", "skyblue1", "indianred1", "mediumpurple1", "olivedrab3", "pink", "#FFED6F", "mediumorchid3", "green" , "#9999CC", "#663366", "#999966", "#9999FF", "seashell1", "skyblue1", "yellow", "red", "olivedrab3", "salmon", "#FFED6F", "mediumorchid3", "gray50", "tan1",  "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "deeppink4","wheat1", "#66CCCC", "forestgreen", "yellow4", "darkorange3"))(35)
barplot(y3, space=0, cex.names=0.15, col=pal, las=2, legend.text = TRUE, lwd = 0.25,
        args.legend = list(x = "topright", bty = "n", inset=c(-0.005, -0.05), cex=0.15))

#Save the image as a pdf to edit in Adobe illustrator

#-----------------------------------------------
#if you want to make in ggplot for prettier colours
#must convert data to long form for ggplot with tidyr
#first transpose so samples are rows, genera by columns
y3t<-t(y3)
#make the sample ID (rownames) as column 1 in the data frame
y3t2 <- data.frame(sampleID = row.names(y3t), y3t)

#convert to long format
library(tidyr)
y3.long2 <- pivot_longer(y3t2, cols=-1, names_to = "Genus", values_to = "Abundance")

x<-y3.long2
#lock in the (factor) order of y3.long2 otherwise ggplot will just plot in alphabetical order
x$sampleID <- factor(x$sampleID, levels = unique(x$sampleID))

ggplot(x, aes(fill=Genus, y=Abundance, x=sampleID)) + 
    geom_bar(position="stack", stat="identity")
