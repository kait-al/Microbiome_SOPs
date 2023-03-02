# Outlier and contaminant detection in microbiota data

These are example R scripts to determine outlier samples and potential contaminant SVs.

---
## Determining outlier samples with CoDaSeq

First, read in your counts table. Here we are working with "filtered_counts" which is generated with [04_filter_and_pca.R](04_filter_and_pca.R). Load in the required libraries:
```r
library(zCompositions)
library(CoDaSeq)

# List samples based on cohorts of interest
# Stone-formers Pre-Op urine samples
S<-c("X101_PU","X102_PU","X103_PU","X104_PU","X105_PU","X106_PU","X107_PU","X108_PU","X110_PU","X111_PU","X112_PU","X113_PU","X114_PU","X115_PU","X116_PU","X117_PU","X118_PU","X119_PU","X120_PU","X121_PU","X122_PU","X123_PU","X124_PU","X125_PU","X126_PU","X127_PU","X128_PU","X129_PU","X130_PU","X131_PU","X132_PU","X133_PU","X134_PU","X135_PU","X136_PU","X137_PU","X138_PU","X139_PU","X140_PU","X141_PU","X142_PU","X143_PU","X144_PU","X145_PU","X146_PU","X147_PU","X148_PU","X149_PU","X150_PU","X151_PU","X152_PU","X153_PU","X154_PU","X155_PU","X156_PU","X157_PU","X158_PU","X159_PU","X160_PU","X161_PU","X162_PU","X163_PU","X164_PU","X165_PU","X166_PU","X167_PU","X168_PU","X169_PU","X170_PU","X171_PU","X172_PU","X173_PU","X174_PU","X175_PU","X176_PU","X177_PU","X178_PU","X179_PU","X180_PU","X181_PU","X182_PU","X183_PU","X184_PU")

# Healthy control urine samples
H<-c("X450_PU","X451_PU","X452_PU","X453_PU","X454_PU","X455_PU","X457_PU","X458_PU","X459_PU","X460_PU","X462_PU","X463_PU","X464_PU","X465_PU","X467_PU","X469_PU","X470_PU","X471_PU","X472_PU","X474_PU","X475_PU","X476_PU","X477_PU","X478_PU","X479_PU")
```

To CLR-transform our data, we can't have zero. We will impute zeroes with the Bayesian-Multiplicative replacement from the zCompositions package. Then we will clr-transform from the CoDaSeq package.
```r
# filtered_counts has samples as rows
czm <- cmultRepl(filtered_counts, label = 0, method = "CZM", output="p-counts")
# czm has samples as rows
d.n0.clr <- codaSeq.clr(czm, samples.by.row=TRUE)
#d.n0.clr has samples as rows
```
We will look at outliers in our cohort groups, as well as in the dataset as a whole. Look at the histogram of the proportional variance for each output.
```r
# get proportional variance per sample in each group
group.H <- d.n0.clr[c(H),]
group.S <- d.n0.clr[c(S),]

# Healthy controls
pvar.H <- codaSeq.outlier(group.H)
pvar.H # look at the entire output
pvar.H$bad # look at just the outliers
# "X454_PU"

# Stone formers
pvar.S <- codaSeq.outlier(group.S)
pvar.S$bad

# Whole dataset
pvar.all <- codaSeq.outlier(d.n0.clr)
```
So there are several potential outliers, the most significant of which are X139_PU and X454_PU. We have statistical substantiation to remove these samples.

---
## Determining potential contaminating features with decontam. See the decontam vignette [here](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) for more info.

This utilizes phyloseq. Generate a phyloseq object with your counts, taxa, and metadata table as was done previously in [05_phyloseq_alpha_diversity.R](05_phyloseq_alpha_diversity.R). This physeq object is based on the unfiltered raw counts. 

The decontam R package has multiple modes of contaminant detection, based on whether the sample is a true sample or a negative control, and based on the input DNA concentration. Your metadata table must contain a column called 'quant_reading' where you've included the input DNA concentration (often provided in an excel spreadsheet from LRGC). It must also contain a column called 'Sample_or_Control', where you specify your true samples as "Sample" and your negative controls (PCR blanks and DNA extraction blanks) as "Control". The negative controls were removed from the dataset in [example_files](example_files) for simplicity, but scripts are included below that you can use with your own dataset that includes your negative controls.

Load the required libraries. Then we will inspect the library sizes:
```r
library(phyloseq)
library(dplyr)
library(decontam)
library(ggplot2)

#physeq was generated previously in '05_phyloseq_alpha_diversity.R'
physeq
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 935 taxa and 108 samples ]
#sample_data() Sample Data:       [ 108 samples by 60 sample variables ]
#tax_table()   Taxonomy Table:    [ 935 taxa by 6 taxonomic ranks ]

df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data frame
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
```
Again, this dataset contains only true samples. The library sizes range from 1016 to 455151.

You can use the following scripts on your own dataset (that includes controls) to determine potential contaminants based on Prevalence mode (whether SVs are found in true samples or negative controls).

```r
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_Control == "Control"

contamdf.prev <- isContaminant(physeq, conc = NULL, neg = "is.neg", method = "prevalence")
table(contamdf.prev$contaminant)

which(contamdf.prev$contaminant)

contamdf.prev01 <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
# you can change the threshold to 0.5 for a more aggressive classification of contaminants


ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev01$contaminant)

# plot this df
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# write a table of the output
write.table(contamdf.prev01, file="contam.prev.01.txt", sep="\t",  quote=F) 
```

Now we will determine potential contaminants with on the Frequency mode (based on DNA concentration). If you don't set a threshold, the default is p < 0.1.
```r
contamdf.freq <- isContaminant(physeq, method="frequency", conc="quant_reading", threshold = 0.05)
head(contamdf.freq)
# the output is a dataframe with columns including 'p' which is the probability for classifying contaminants

table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

#make a table to look at the contaminants in order of the lowest p 
order <- contamdf.freq[order(contamdf.freq$p, rev(contamdf.freq$p), decreasing = FALSE),]
#look at the first 20 rows
order[1:20,]
```
So 7 of the 935 SVs are considered contaminants with this mode. Apparently SV_96 is the most likely to be an outlier, so we will look at it on a plot beside SV_0 which is unlikely to be an outlier.
```r
plot_frequency(physeq, taxa_names(physeq)[c(1,97)], conc="quant_reading") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")
```
![alt text](https://github.com/kait-al/Microbiome_SOPs/blob/main/Amplicon_SOPs/images/decontam.freq.jpg)

The dashed black line shows the model of a noncontaminant sequence feature for which frequency is expected to be independent of the input DNA concentration. The red line shows the model of a contaminant sequence feature, for which frequency is expected to be inversely proportional to input DNA concentration, as contaminating DNA will make up a larger fraction of the total DNA in samples with very little total DNA.
```r
# Make data.frame of prevalence in positive and negative samples, instead colouring by freq this time instead of prev
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.freq$contaminant)

# plot this df, how do the contaminants compare to those detected with prevalence mode?
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
```

You can also determine potential contaminants with the method "combined", which considers both prevalence and frequency with Fisher's method.
```r
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_Control == "Control"
contamdf.comb <- isContaminant(physeq, method="combined", conc="quant_reading", neg="is.neg")

# Make data.frame of prevalence in positive and negative samples, instead colouring by combined this time 
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.comb$contaminant)

# plot this df, how do the contaminants compare to those detected with prevalence mode?
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
```