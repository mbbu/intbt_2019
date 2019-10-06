#SCRIPT PURPOSE: DEMONSTRATION OF 16S DOWNSTREAM ANALYSES (USING THE 16S PACKAGES metagenomeSeq, phyloseq, vegan)
#AUTHOR: KATIE LENNARD

# Load custom functions, set working directories, import data --------------------------
source("microbiome_custom_functions_tutorial.R")
#Read in ASV counts
ASVs <- readRDS("ASVs_counts.RDS")#Default options with primer trimming & new refseq_rdp DB
str(ASVs)
dim(ASVs)
#Read in taxonomy table
taxa <- readRDS("ASVs_taxonomy.RDS")
dim(taxa)
#
identical(rownames(ASVs),rownames(taxa))#TRUE
#Assign user-friendly ASV IDs to replace sequences
head(rownames(ASVs))
seqs <- rownames(ASVs)
ASV.IDs <- paste0("ASV",c(1:length(seqs)))
#Named vector:
names(seqs) <- ASV.IDs
head(seqs)
seq_lens <- nchar(seqs)
seq_lens
plot(density(seq_lens))
#Merge ASV table and taxonomic table as a phyloseq object
phy <- phyloseq(otu_table(ASVs,taxa_are_rows=TRUE),tax_table(taxa))
identical(taxa_names(phy),rownames(ASVs))#TRUE
taxa_names(phy) <- names(seqs)
str(phy)
tax_df <- data.frame(tax_table(phy))
tax_df$seq_len <- seq_lens
head(tax_df)
dim(tax_df)
length(which(is.na(tax_df[,"Species"])))# if GG ~ 32% assigned | If RefSeq-RDP ~ 70% assigned at species-level (3048/4415)

# Import metadata --------------------------------------------------
meta <-  read.table("practice.dataset1.metadata.tsv", sep = "\t", header =TRUE, row.names=1)
head(meta)
##       Dog Treatment
## Dog1    B         2
## Dog2    G         3
## Dog3    K         3
## Dog8    B         4
## Dog9    G         0
## Dog10   K         4

rownames(meta)

##  [1] "Dog1"  "Dog2"  "Dog3"  "Dog8"  "Dog9"  "Dog10" "Dog15" "Dog16"
##  [9] "Dog17" "Dog22" "Dog23" "Dog24" "Dog29" "Dog30" "Dog31"

head(sample_names(phy))

## [1] "Dog10" "Dog15" "Dog16" "Dog17" "Dog1"  "Dog22"

length(sample_names(phy))#15

## [1] 15

length(rownames(meta))#15 (check if same number of samples in .biom file and metadatafile)

## [1] 15

length(intersect(rownames(meta),sample_names(phy)))#15 (check that the sample names match in all cases)

## [1] 15

#--------------------
# Assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
sample_data(phy) <- meta 
identical(sample_names(phy),colnames(otu_table(phy)))#see
#--------------------
str(sample_data(phy)) # Need to change treatment column to factor variable

# 'data.frame':	15 obs. of  2 variables:
#   Formal class 'sample_data' [package "phyloseq"] with 4 slots
# ..@ .Data    :List of 2
# .. ..$ : Factor w/ 3 levels "B","G","K": 1 3 1 2 3 2 1 2 3 1 ...
# .. ..$ : int  2 4 1 4 0 3 3 1 2 0 ...
# ..@ names    : chr  "Dog" "Treatment"
# ..@ row.names: chr  "Dog1" "Dog10" "Dog15" "Dog16" ...
# ..@ .S3Class : chr "data.frame"

sample_data(phy)[,"Treatment"] <- as.factor(unlist(sample_data(phy)[,"Treatment"]))
saveRDS(phy, "dada2_dog_stool.RDS") 
readRDS("dada2_dog_stool.RDS")
# Exploratory analysis ----------------------------------------------------

#Let's explore the data
reads <- sample_sums(phy) #number of reads per sample
reads
#Set (arbitrary) cutoff for number of acceptable reads/sample. Note this depends on your sequencing platform (MiSeq vs. HiSeq) and sample type
length(which(reads<5000)) #all samples have at least 5000 reads
#Let's standardize sample read count so we can compare between samples:
total = median(sample_sums(phy))#find median sample read count
standf = function(x, t=total) round(t * (x / sum(x)))#function to standardize to median sample read count
M.std = transform_sample_counts(phy, standf)#apply to phyloseq object
ntaxa(M.std)#4415
sample_sums(M.std)
#----
#Lets plot alpha (within sample) diversity with the plot_richness function from the phyloseq package
p <- plot_richness(M.std,x = "Dog",color = "Treatment",measures=c("Shannon"), 
				title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=16,face="bold"))+geom_point(size=5)
p

#Note: Regarding the warning on the absence of singletons, please see https://github.com/benjjneb/dada2/issues/214 (and ignore warning)

pdf("alpha_diversity_by_dog_treatment.pdf")
p
dev.off()

est <- estimate_richness(M.std, split = TRUE, measures = c("Shannon"))
temp <- cbind(est,sample_data(M.std)[,"Dog"])
head(temp)
# Shannon Dog
# Dog1    5.638   B
# Dog10   4.944   K
# Dog15   5.172   B
# Dog16   5.793   G
# Dog17   5.402   K
# Dog2    5.590   G
aggregate(temp$Shannon, list(temp$Dog), mean)

t <- kruskal.test(temp[,1]~temp[,2])
t
# Kruskal-Wallis rank sum test
# 
# data:  temp[, 1] by temp[, 2]
# Kruskal-Wallis chi-squared = 4.6, df = 2, p-value = 0.1
#OR, instead of plotting alpha diversity by dog we can plot by treatment

p <- plot_richness(M.std,x = "Treatment",color = "Dog",measures=c("Shannon"), 
				title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=16,face="bold"))+geom_point(size=5)
p

pdf("alpha_diversity_by_treatment_dog.pdf")
p
dev.off()
temp <- cbind(est,sample_data(M.std)[,"Treatment"])
head(temp)
# Shannon Treatment
# Dog1    5.638         2
# Dog10   4.944         4
# Dog15   5.172         1
# Dog16   5.793         4
# Dog17   5.402         0
# Dog2    5.590         3

t <- kruskal.test(temp[,1]~temp[,2])
t
# Kruskal-Wallis rank sum test
# 
# data:  temp[, 1] by temp[, 2]
# Kruskal-Wallis chi-squared = 2.8, df = 4, p-value = 0.6


# Beta diversity ----------------------------------------------------------
#Now that we've looked at alpha diversity of the unfiltered dataset, lets look at beta diversity on the merged taxonomies - WHY? Because e.g. Bray-Curtis dissimilarity
#sees ASV1 and ASV2 as two entirely different entities, even though they acdtually belong to the same species. We want to summarize biologically relevant data.
#Let's merge our ASVs at their lowest available taxonomic annotation level. Why? Look at tax_df - can you see there are several ASVs for a given species? Strain level variation. 
#For our purpose we would like to examine at species-level or higher. We therefore first merge and THEN filter

#It's useful to collapse (sum) ASVs at their lowest available taxonomic annotations:
#For this we use the tax_glom.kv() function from the 'microbiome_custom_functions_tutorial.R' script loaded at the beginning of this script
M.phy <- tax_glom.kv(M.std)

## [1] "Removing phylogenetic tree"
## [1] "There are now xx merged taxa"
ntaxa(M.phy)
## [1] 124

#Let's filter out taxa that are only present at very low numbers in a small minority of samples
#The filter below retains only OTUs that are present at at least 10 counts at least 20% of samples OR that have a total relative abundance of at least 0.1% of the total number of reads/sample
M.f = filter_taxa(M.phy,function(x) sum(x > 10) > (0.2*length(x)) | sum(x) > 0.001*total, TRUE)
ntaxa(M.f)#106

#Now, let's look at beta (between-sample) diversity, using non-metric multidimensional scaling (NMDS)
?ordinate #from the phyloseq package performs ordination of phyloseq data
set.seed(20) #select any start seed, in this case 2 (ensures reproducibility)
ord.BC <- ordinate(M.f, "NMDS", "bray", k=2, trymax=1000) # stress=0.07
stressplot(ord.BC)
goodness(ord.BC)
ord.BC
color = c("Treatment")
shape = c("Dog")
title=c("NMDS of 16S microbiome,Bray-Curtis distance,k=2")
MDS = plot_ordination(M.f, ord.BC, color = color,shape=shape, 
		title = title)
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
		theme_bw()+labs(color=color, shape=shape)+geom_point(size=5)

MDS.1

pdf("NMDS_Dogs_treatment_Bray_Curtis.pdf",6,5)
MDS.1
dev.off()

#NB what happens if we use unmerged ASVs instead of merged ASVs to calculate Beta diversity?

test.BC <- ordinate(M.std, "NMDS", "bray", k=2, trymax=1000) # stress=0.07
# Run 1000 stress 2.170209e-06 
# ... Procrustes: rmse 0.2549385  max resid 0.3839214 
# *** No convergence -- monoMDS stopping criteria:
#   9: no. of iterations >= maxit
# 991: stress < smin
# Warning messages:
#   1: In metaMDS(veganifyOTU(physeq), distance, ...) :
#   stress is (nearly) zero: you may have insufficient data
# 2: In postMDS(out$points, dis, plot = max(0, plot - 1), ...) :
#   skipping half-change scaling: too few points below threshold
stressplot(test.BC)
goodness(test.BC)

# Heatmap -----------------------------------------------------------------
filename <- c("heatmap_merged_taxa")

main <- paste("Merged taxa, Bray-Curtis distance")

f = paste0(filename,".pdf")
# Color specification for column annotations above heatmap:
D.cols = c("B"="#CC79A7","G"="#56B4E9","K"="#F0E442")
colours = list(Dog=D.cols)

# Create sample pairwise distance matrix (Bray-Curtis distance) and cluster with hclust():
set.seed(2)
diss <- distance(M.f,method = "bray", type = "samples")

clust.res<-hclust(diss)
sample.order = clust.res$order
# Heatmap is output to file (the heatmap.k function can be found in the 'microbiome_custom_functions.R' script)
hm = heatmap.k(physeq= M.f, annot.cols = c(1,2), main = main,filename = f,colours=colours,cexRow = 0.5,Colv = sample.order,labrow = TRUE)  


# Barplot -----------------------------------------------------------------

#Let's also explore the data at Genus-level as barplot by Dog
#NOTE: We use the unmerged data as input here as we now wish to merge at Genus level and not the lowest available taxonomic level
level = "Genus"
count = 500
perc = 0.25
# Barplot will be written to file (the bar.plots function can be found in the 'microbiome_custom_functions.R' script)
barplot = bar.plots(physeq = M.std,cat = "Dog",level = level, count = count, perc = perc, outDir=outDir, 
		filen = 'Barplots_by_Dog')
print(barplot)

# Differential abundance testing ------------------------------------------
#For diferential abundance testing we use raw reads (not standardized) as input because metagenomeSeq has its 
#own internal scaling factor based on total read counts/sample

#So lets merge the raw reads at lowest available taxonomy - again we're interested in biologically interpretable units, not strain-level variation
Mraw = tax_glom.kv(phy)
ntaxa(Mraw)#124
Mraw.f = filter_taxa(Mraw,function(x) sum(x > 10) > (0.2*length(x)) | sum(x) > 0.001*total, TRUE)
ntaxa(Mraw.f)#107

#Let's subset the metagenomeSeq object to eliminate samples from dog K (we would like to compare dogs G & B)
sub.index <- sample_names(Mraw.f)[sample_data(Mraw.f)[,"Dog"] != "K"]
phy.temp <- prune_samples(sub.index, Mraw.f)
nsamples(phy.temp)
## [1] 10
ntaxa(phy.temp)
#107
#LEt's check if there are ASVs where our subset of samples now have only zero counts (i.e. ASVs that are exclusively found in dog K)
length(which(rowSums(otu_table(phy.temp))==0))#1
keep <- which(rowSums(otu_table(phy.temp))!=0)
#LEts remove these empty rows
phy.temp <- prune_taxa(names(keep),phy.temp)
ntaxa(phy.temp)#106

#The function we'll use to do differential abundance testing was built around metagenomeSeq's fitZig() function
#to build the model and MRfulltable() to extract the results. This customized function also filters the results by fold-change, percent presence across samples and adjusted p-values
#Results are summarized as a heatmap and .csv table

#There are several parameters that we can specify in the super.fitZig.kv() function including:
#ARGUMENTS:
#physeq = phyloseq object with raw counts but pruned by read count (e.g exclude samples with < 10 000 reads)
#factor = the factor of interest (This should be a column name in sample_data(physeq))
#FileName = specify comparison made, not full file path (e.g. "M.prune_PICRUSt_groups_C_vs_A")
#p = adjusted p-value threshold to use when filtering significant results.
#default p=0.05
#perc = threshold for filtering by fraction of +ve ssamples for a given OTU in either group (e.g perc = 0.3 will keep only OTUs where at least one of the two groups have ? 30% of samples +ve for that OTU) 
#default perc = 0.5
#FC = absolute coefficient threshold to use for filtering (default = 1.5)
#main = title of heatmap
#subt = subtitle for heatmap (most commonly filtering settings e.g."FDR < 0.05,|coeff| >= 1.5, >30%+ in either group"
#heatmap.describer = heatmap filename segment specifying clustering or not as well as type of annotation (OTU ID/taxonomy) (e.g."category_ordered_tax_annot.pdf" )
#order = should heatmap columns be sorted manually? TRUE/FALSE, default=TRUE	
#rows.sortby = would you like to sort the resulting heatmap by a column from the rsults table e.g. "adjPavalues" or "coeff"
#labrow = should the heatmap be annotated with taxonomic annotations for OTU ids (TRUE/FALSE; TRUE = taxonomic annotations; FALSE = OTU IDs)
#extra.cols = any extra columns that you would like to plot on the heatmap e.g. c("feeding","HIV_exposure")
#cexCol=scaling factor for heatmap column labels (usually sample names)
#cexRow=scaling factor for heatmap row labels (usually taxon names)


#Note that we do not HAVE to specify each and every one of these parameters. Some have default values that will be used if we don't specify anything.
b = super.fitZig.kv(physeq = phy.temp,factor = "Dog",outDir = outDir,FileName =c("1_25FC_0.2_Dog_GvsB_ASVs"),
		heatmap.descriptor=c("tax_annot"), main=c("Dog G vs. B, merged taxa"), subt=c("subt = FDR < 0.05,|coeff| >= 1.25, >20%+ in either group"), 
		ordered=TRUE, p=0.05, FC = 1.25, perc=0.2, extra.cols = c("Treatment"))

b

#///////////////////////////////
