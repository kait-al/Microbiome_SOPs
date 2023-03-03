#Aldex ttest of differential abundance between 2 groups

library(ALDEx2) 

counts <- read.table("data/cutadapt_counts.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)
#or use filtered_counts

#transpose
tc<-t(counts)

#specify your groups

#Stone former samples begin with X1
SF<-grep("X1", colnames(tc), value = TRUE)
#Healthy control samples begin with X4
HC<-grep("X4", colnames(tc), value = TRUE)

#if groups can't be grep'ed out by sample ID, provide a list of column names for each group
#for example:
#SF<-c("X101_PU","X102_PU","X103_PU","X104_PU","X105_PU","X106_PU","X107_PU","X108_PU","X109_PU","X110_PU","X111_PU","X112_PU","X113_PU","X114_PU","X115_PU","X116_PU","X117_PU","X118_PU","X119_PU","X120_PU","X121_PU","X122_PU","X123_PU","X124_PU","X125_PU","X126_PU","X127_PU","X128_PU","X129_PU","X130_PU","X131_PU","X132_PU","X133_PU","X134_PU","X135_PU","X136_PU","X137_PU","X138_PU","X139_PU","X140_PU","X141_PU","X142_PU","X143_PU","X144_PU","X145_PU","X146_PU","X147_PU","X148_PU","X149_PU","X150_PU","X151_PU","X152_PU","X153_PU","X154_PU","X155_PU","X156_PU","X157_PU","X158_PU","X159_PU","X160_PU","X161_PU","X162_PU","X163_PU","X164_PU","X165_PU","X166_PU","X167_PU","X168_PU","X169_PU","X170_PU","X171_PU","X172_PU","X173_PU","X174_PU","X175_PU","X176_PU","X177_PU","X178_PU","X179_PU","X180_PU","X181_PU","X182_PU","X183_PU","X184_PU")

aldex.in<-d[,c(HC, SF)]

conds<-c(rep("HC", length(HC)), rep("SF", length(SF)))

#get the clr values
#this is the main ALDEx function for all downstream analyses
#mc.samples=128 is often sufficient
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)

#perform t-test (both Welches and Wilcoxon, plus a Benjamini-Hochberg multiple test correction)
x.tt <- aldex.ttest(x, paired.test=FALSE)

#estimate effect size and the within and between condition values
#include indiv. samples or not
x.effect <- aldex.effect(x)

#merge the data
x.all <- data.frame(x.tt, x.effect)

#get features passing significance
sig <- which(x.tt$we.eBH < 0.05)
sig

#significant BH Wilcoxon > 0.1 (wi.eBH)
#several have an effect size > |0.5| (significant!!)

#write a .txt with the results
write.table(x.all, file="aldex_HCvSF_output.txt", sep="\t", quote=F, col.names=NA)
