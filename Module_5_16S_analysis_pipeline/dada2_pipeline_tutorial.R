# ================================================================================
# 16S rRNA pipeline using DADA2
#
# Tutorial by Imane Allali
#
# Available online at: https://iallali.github.io/DADA2_pipeline/16SrRNA_DADA2_pipeline.html
# ================================================================================

# load package
library(dada2); packageVersion("dada2")

# set path to data
MY_HOME <- Sys.getenv("HOME")
data <- paste(MY_HOME, "/dada2_tutorial_dog/dog_samples", sep='')  # change the path
list.files(data)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
dataF <- sort(list.files(data, pattern="_R1.fastq", full.names = TRUE))
dataR <- sort(list.files(data, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
list.sample.names <- sapply(strsplit(basename(dataF), "_"), `[`, 1)
list.sample.names

####################################################################################
# STEP 1: Quality Control  on the Raw Data
####################################################################################

# visualize the quality profiles of the datasets

# forward
plotQualityProfile(dataF[1:3])

# reverse
plotQualityProfile(dataR[1:3])

####################################################################################
# STEP 2: Filter and trim the raw data
####################################################################################

# Place filtered files in filtered/ subdirectory
filt.dataF <- file.path(data, "filtered", paste0(list.sample.names, "_F_filt.fastq.gz"))
filt.dataR <- file.path(data, "filtered", paste0(list.sample.names, "_R_filt.fastq.gz"))
names(filt.dataF) <- list.sample.names
names(filt.dataR) <- list.sample.names

out <- filterAndTrim(dataF, filt.dataF, dataR, filt.dataR, truncLen=c(290,275),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

# visulize
head(out)

####################################################################################
# STEP 3: Learn the error rates
####################################################################################
# DADA2 pipeline uses the learnErrors method to learn the error model from the data

#forward
errF <- learnErrors(filt.dataF, multithread=TRUE)

# reverse
errR <- learnErrors(filt.dataR, multithread=TRUE)

# visualize the estimated error rates for the forward reads
plotErrors(errF, nominalQ=TRUE)

# visualize the estimated error rates for the reverse reads
plotErrors(errR, nominalQ=TRUE)

####################################################################################
# STEP 4: Sample Inference
####################################################################################

# run dada2 algorithm on trimmed and filtered forward and reverse reads

# forward reads
dadaF <- dada(filt.dataF, err=errF, multithread=TRUE)

#reverse reads
dadaR <- dada(filt.dataR, err=errR, multithread=TRUE)

# visualize the dada2-class object
dadaF[[1]]

####################################################################################
# STEP 5: Merge the Paired reads
####################################################################################

merge.reads <- mergePairs(dadaF, filt.dataF, dadaR, filt.dataR, verbose=TRUE)

# inspect the merger data.frame from the first sample
head(merge.reads[[1]])

####################################################################################
# STEP 6: Construct the Sequence table (amplicon sequence variant table (ASV))
####################################################################################

seqtab <- makeSequenceTable(merge.reads)
dim(seqtab)

# inspect the distribution of sequence lengths
table(nchar(getSeqences(seqtab)))

####################################################################################
# STEP 7: Remove Chimeras
####################################################################################
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# check the dimensions
dim(seqtab.nochim)

# percentage of sequences that are not chimeras
sum(seqtab.nochim)/sum(seqtab)

####################################################################################
# STEP 8: Track Reads through the DADA2 Pipeline
####################################################################################

getN <- function(x) sum(getUniques(x))
track.nbr.reads <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(merge.reads, getN), rowSums(seqtab.nochim))

colnames(track.nbr.reads) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track.nbr.reads) <- list.sample.names
head(track.nbr.reads)

####################################################################################
# STEP 9: Assign Taxonomy
####################################################################################

# Download RefSeq-RDP database from https://zenodo.org/record/3266798/files/RefSeq-RDP16S_v3_May2018.fa.gz
# assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, paste(MY_HOME, "/dada2_tutorial_dog/RefSeq-RDP16S_v3_May2018.fa.gz", sep=''), multithread=TRUE) # change the path

# visulise taxonomy after assignment
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# save the ASVs table in your working directory.
write.csv(taxa, file="ASVs_taxonomy.csv")
saveRDS(taxa, "ASVs_taxonomy.rds")

# save the ASV table
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
count.asv.tab <- t(seqtab.nochim)
row.names(count.asv.tab) <- sub(">", "", asv_headers)
write.csv(count.asv.tab, file="ASVs_counts.csv")
saveRDS(count.asv.tab, file="ASVs_counts.rds")

####################################################################################
# STEP 10. Alignment
####################################################################################

# Before constructing the phylogenetic tree, we should do an anlignment
# uses the DECIPHER package
library(DECIPHER)

# We use the seqtab (sequence table as an input) that we got from step number 6
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

####################################################################################
# STEP 11: Construct Phylogenetic Tree
####################################################################################

# We use the phangorn package
library(phangorn)

# first construct a neighbor-joining tree then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree as a starting point
phang.align <- phyDat(as(alignment, "matrix"), type="DNA") 
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

# save the fitGTR file
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR, "phangorn.tree.RDS")








