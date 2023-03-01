# Pruning and filtering your counts tables
#### This was adapted from [this very helpful tutorial](https://github.com/jmmack/16S/blob/master/manipulating_counts_table.md) by Jean Macklaim.

Here are some examples on how to filter and inspect your data. Inspect your dataset often!! Make sure it looks right after you manipulate it and before you proceed.

---
## Inspection

Import your raw, unfiltered counts table into R. The counts table should have samples as columns and features as rows. Features can be amplicon sequence variants (SVs), genes, functional pathways, other taxonomic levels, etc. This SV table can be found in [example_files](example_files). 
```r
d<-read.table("dada2_nochim_tax.txt", sep="\t", quote="", header=T, row.names=1)
# include comment.char="" if your header line starts with #
```

**Only** if your table has samples as rows, transpose it to use the rest of the code as-is:
```r
d <- t(d)
```

Check your table dimensions (number of rows and columns):
```r
dim(d)
#[1] 935 284
```

Look at your column and row names:
```r
colnames(d)
rownames(d)
```

Look at the first 5 rows and columns of d:
```r
d[1:5,1:5]
#      X000_CYSTINE X101_OU X101_PU X101_S_L X101_S_R
# SV_0        10441   11004    2732       12      852
# SV_1            0    1038   92126      359     5000
# SV_2          182       0     262      193        0
# SV_3          710       0       0        0        0
# SV_4           25     116     121        0     2028

```

Check the sparsity of your dataset (number of zeroes):
```r
sum(d == 0)
# 252083

sum(d != 0)
#13457
```

Like most microbiota data, this is a very sparse dataset! There are many SVs that are detected in only a handful of samples. 

Overall, we have 935 SVs and 283 samples. The last column (284) is the taxonomy string.

If, like this example, your taxonomy string is in the last column (here the column name is "tax.vector"), we have to remove it before proceeding. So that we dont lose it altogether, we will temporarily move it to another object called tax:
```r
tax<-d$tax.vector
#get only the count columns
dm<-d[,1:ncol(d)-1]
dim(dm)
#935 283
```
---

## Filter samples (columns)

Keep only columns with >1000 reads:
```r
i <- (colSums(dm) <=1000)
d.s <- dm[, !i]
dim(d.s)
# 935 239
ncol(dm)-ncol(d.s)
# 44
```

So we went from 283 to 239 samples (44 removed).

---

## Filter SVs (rows) based on different *frequency* cutoffs.

First calculate frequency:
```r
d.freq <- apply(d.s, 2, function(x){x/sum(x)})
```

Think about your own dataset here. Based on your biological hypothesis, do you care more about rare or abundant taxa? Are samples derived from high (gut, vaginal, oral) or low (urine, tissue) bacterial abundance communities? How much of a role will contamination play? You can experiment with different cutoffs.
I commonly use 0.01, 0.05, and 0.001. 

Here we will start by keeping SVs with a frequency of > 0.01 (1%) in *any* sample: 
```r
d.0 <- d.s[apply(d.freq, 1, max)>0.01,]
dim(d.0)
# 205 239 
#So we went from 935 SVs to 205.
```

Now we will keep SVs with a frequency of > 0.01 (1%) in *every* sample:
```r
d.1 <- d.s[apply(d.freq, 1, min)>0.01,]
dim(d.1)
# 0 239
```

While that last filtering technique might make sense in gut or vaginal samples, these are urinary samples, and there is *not a single SV* present in at least 1% abundance in *every* sample.

---

## Filter SVs (rows) based on a *read count* cutoff

Keep only SVs that have a *total* read count (row sum) across all samples of 1000.
If you have a lot of samples (> 500 columns), you could increase this to 10000 reads. For a dataset with fewer samples, you could use 100.
```r
count = 1000
d.2 <- data.frame(d.s[which(apply(d.s, 1, function(x){sum(x)}) > count),])
dim(d.2)
# 211 239
# So we went from 935 to 211 SVs
```

Keep only SVs that have a *mean* count of at least 1 across all samples:
```r
count = 1
d.3 <- data.frame(d.s[which(apply(d.s, 1, function(x){mean(x)}) > count),])
dim(d.3)
# 352 239
# So we went from 935 to 352 SVs
```

Discard SVs if it is a zero in half or more of the samples:
```r
cutoff = .5
d.4 <- data.frame(d.s[which(apply(d.s, 1, function(x){length(which(x != 0))/length(x)}) > cutoff),]) 
dim(d.4)
# 23 239
# So we went from 935 to 23 SVs (97% of SVs are 0 in more than half of the samples!)
```

---

## Select samples (columns) by name

Make a vector of the column names to keep for downstream analysis. These can be separated into comparison groups of interest.
```r
# Stone-formers Pre-Op urine
# n = 84, only 5 shown here for simplicity
SPU<-c("X101_PU","X102_PU","X103_PU","X104_PU","X105_PU")

# Stone-formers OR urine
# n = 59 after filtering, only 5 shown here for simplicity
SOU<-c("X101_OU","X102_OU","X103_OU","X104_OU","X105_OU")

# Healthy control urine samples
# n = 25 after filtering, only 5 shown here for simplicity
HCU<-c("X450_PU","X451_PU","X452_PU","X453_PU","X454_PU")

# If you want a separate dataframe of just these samples, concatenate them together. Use whichever filtered table you want (d, d.s, d.0 to d.4).
df1 <- d.3[,c(SPU, SOU, HCU)]
dim(df1)
# 352 15
```

Make a vector of the column names to exclude for downstream analysis. For example, remove your negative and positive controls if they are significantly distinct from your clinical samples (distinct on a biplot, or statistically distinct for example with envfit/permanova).
```r
remove<-c("DNAB1","DNAB3","PCRB1","PCRB2","PCRB3","SPIKE1_1","SPIKE1_2","SPIKE1_3","SPIKE2_1","SPIKE2_2","SPIKE2_3")
length(remove)
#11
df2<-d.3[, !names(d.3) %in% remove]
dim(df2)
#352 228
```

Grep to include/exclude samples based on pattern matching. All the Healthy control urine samples start with "X4", whereas samples starting with "X1" are from stone formers.
```r
# make a new data frame of all stone former samples
df.sf <- d.3[, grep("X1", colnames(d.3))]
dim(df.sf)
# 352 195

# make a new data frame of all healthy control samples
df.hc <- d.3[, grep("X4", colnames(d.3))]
dim(df.hc)
# 352 25

# make a new dataframe of all samples except healthy control samples 
removeHC<-grep("X4", colnames(d.3)) 
df.no.hc <- d.3[-c(removeHC)]
dim(df.no.hc)
# 352 214
# This df still includes positive and negative controls which is why dims are more than df.sf
```
---

Although these filtering steps were all separately performed on d.freq and d.s, you could compound any of these filtering steps together. 

***Important***: As soon as you remove samples (or select a subset of samples), you are removing SV read counts. Therefore, the previous filters will no longer apply and you should re-filter. A 1% abundance filter from all samples is not the same as a 1% abundance filter on a subset of samples because the counts are different (the sample subset might not contain any of SV#X even though SV#X passed the 1% filter in the original dataset)!

---
## Final steps

When you're done, if you want to put your taxonomy info back in:
```r
d.3$tax.vector = d$tax.vector[match(rownames(d.3), rownames(d))]
```
Or, you could make taxonomy your rownames earlier on, and remove the tax column. This only works without duplicate taxa:
```r
rownames(d)<-d$tax.vector
d$tax.vector <- NULL
```
Or, you could make taxonomy PLUS OTU number your rownames:
```r
rownames(d)<-paste(rownames(d), d$tax.vector, sep=";")
d$tax.vector <- NULL
```

Write your final table with a descriptive name to use in your downstream analyses.
```r
write.table(d, file="filtSV_01abund_1000minrc.txt", sep="\t", quote=F)
```
