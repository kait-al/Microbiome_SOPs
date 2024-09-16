# PhilR (Phylogenetic Isometric Log-Ratio Transform), an alternative way to get CLR values that are related through phylogenetic tree data

# Load library 
library(ape)
library(castor)
library(philr)

# Load data
counts <- read.table() # Counts table
taxa <- read.table() # Taxa table
meta <- read.table() # Meta data table

# Alignment created using mafft via conda (mafft --auto 16S_sequences.fasta > 16S_sequences_mafft.fasta)
# Tree created using fasttree via conda (FastTree -gamma -gtr -nt 16S_sequences_mafft.fasta > vitk_tree.newick)
tree <- read.tree()

# 16S sequences are extracted from the 8th column of the taxa table
#library(Biostrings)
#### Filter 16S sequences ####
#16S_seqs <- DNAStringSet(taxa[, "Sequence"]) # Make DNAStringSet
#names(16S_seqs) <- rownames(taxa) # Add in the SV names

# Tree from fasttree is not going to be rooted, we will do midpoint rooting via castor
rooted_tree <- root_at_midpoint(tree)

# Perform checks on tree
is.rooted(rooted_tree) # [1] TRUE
is.binary(rooted_tree) # [1] TRUE

# Add internal node labels
rooted_tree <- makeNodeLabel(rooted_tree, method="number", prefix='n')

# You have a few options to proceed here; option 1, czm the counts; option 2, add a pseudocount
# Perform Bayesian-Multiplicative replacement of count zeros
# samples must be as as rows
czm <- cmultRepl(counts, label = 0, method = "CZM", z.warning = 1)


philr_option1 <- philr(as.matrix(czm), 
               rooted_tree,
               part.weights='enorm.x.gm.counts', 
               ilr.weights='blw.sqrt',
               return.all=TRUE)


philr_option2 <- philr(as.matrix(counts), 
                   rooted_tree,
                   part.weights='enorm.x.gm.counts', 
                   ilr.weights='blw.sqrt',
                   pseudocount = 0.1,
                   return.all=TRUE)

# Regardless of the option you pick, the results will be very similar. You will have a new CLR table that facotrs in phylogeny information, which is accessible via philr$x.ilrp
# This can be used to go on and perform PCA (not shown) or PERMANOVA
set.seed(69)
vegan::adonis2(philr$x.ilrp ~ (Group_of_interest), data = meta, method = "euclidean", permutations = 999)

# You can loop through all of your meta data to see if anything is interesting (note this does not take into consideration betadisp)

### PERMANOVA cycle
permanova <- data.frame(Variable = character(), Significance = numeric(), stringsAsFactors = FALSE)
for (col_index in 1:length(colnames(meta))) { 
  col_name <- names(meta)[col_index]
  
  # Construct the formula with the current column
  formula <- as.formula(paste("clr ~", col_name)) #clr should be changed to philr_optionx (you pick)
  
  # Run adonis2 with the current formula
  set.seed(69)
  adonis_result <- adonis2(formula, data = meta, method = "euclidean", permutations = 999, na.action = na.omit)
  
  # Extract only the first column from the 'Pr(>F)' result
  significance_value <- adonis_result$`Pr(>F)`[1]
  
  # Store the results in the data frame
  result_row <- data.frame(
    Variable = col_name,
    Significance = significance_value
  )
  
  permanova <- rbind(permanova, result_row)
}

View(permanova)