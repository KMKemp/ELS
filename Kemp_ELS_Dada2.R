###############################################################
##### Set the directories and load the required libraries #####
############################################################### 
setwd("~/Box/Kemp_MSEW_2020/new_analysis")
base_dir <- '~/Box/Kemp_MSEW_2020/new_analysis/Rcode'
library(dada2)
library(phyloseq)
library(tidyverse)
library(DECIPHER)
library(phangorn)

##########################################################################################
##########################################################################################
##### I. Quality-filter reads and create Amplicon Sequence Variant tables with DADA2 #####
##########################################################################################
##########################################################################################

#####################################################
##### To get started, make the needed variables #####
#####################################################
# One with all the sample names by scanning the "samples" made earlier
samples <- scan("~/Box/Kemp_MSEW_2020/new_analysis/adapter_trimming/samples", what="character")
# Assign path to cutadapt trimmed reads
cutadapt_path <- file.path(base_dir, "adapter_trimming/cutadapt/trimmed_reads") 
# One holding the file names of all the forward reads
forward_reads <- file.path(cutadapt_path, paste0(samples, "_R1_trimmed.fq.gz"))
file.exists(forward_reads) #TRUE
# And one holding the file names of all the reverse reads
reverse_reads <- file.path(cutadapt_path, paste0(samples, "_R2_trimmed.fq.gz"))
file.exists(reverse_reads) #TRUE

#########################
##### Quality check #####
#########################
# Look at the quality of the forward and reverse reads to determine trimming strategy
# https://benjjneb.github.io/dada2/tutorial.html
plotQualityProfile(forward_reads[1:4])
plotQualityProfile(reverse_reads[1:4])

##########################################
##### Perform filtering and trimming #####
##########################################
# Fist make directory a directory called "filtered" within the working directory.
# Then assign a path to place filtered files in the "filtered" subdirectory
filt_path <- file.path(base_dir, "filtered") # Place filtered files in filtered/ subdirectory
# Make variables for forward and reverse filtered reads
filtered_forward_reads <- file.path(filt_path, paste0(samples2, "_R1_filtered.fq.gz"))
filtered_reverse_reads <- file.path(filt_path, paste0(samples2, "_R2_filtered.fq.gz"))
# Run filtering
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(3,3), maxN=0,
                              rm.phix=TRUE, truncQ=2, truncLen=c(220,150), compress=TRUE, multithread=TRUE)
# Make a dataframe to hold filtering stats
filtered_out_df <- data.frame(sample = row.names(filtered_out), filtered_out) %>% mutate(retained = reads.out/reads.in)
print(filtered_out_df)

############################################################
##### Look at the quality after trimming and filtering #####
############################################################
plotQualityProfile(filtered_forward_reads[1:4])
plotQualityProfile(filtered_reverse_reads[1:4])

#################################
##### Learn the error rates #####
#################################
# Forward reads
errF <- learnErrors(filtered_forward_reads, multithread=TRUE)
# Reverse reads
errR <- learnErrors(filtered_reverse_reads, multithread=TRUE)

#####################################
##### Visualize the error rates #####
#####################################
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

############################################################
##### Dereplicate identical reads for faster computing #####
############################################################
derepFs <- derepFastq(filtered_forward_reads, verbose=TRUE)
derepRs <- derepFastq(filtered_reverse_reads, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- samples2
names(derepRs) <- samples2

######################################################
##### Infer the sequence variants in each sample #####
######################################################
# Use the pseudo-pooling option, which makes it easier to resolve rare variants that are present 
# as singletons or doubletone in one sample but are present many times across samples.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")
# Inspect the dada-class object returned by dada:
dadaFs[[1]]
dadaRs[[1]]

##############################
##### Merge paired reads ##### 
##############################
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame for the first sample
print(head(mergers[[1]]))

####################################
##### Construct sequence table #####
####################################
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
print(table(nchar(getSequences(seqtab))))

############################################
##### Plots to examine sequence length #####
############################################
# Plot the sequence variant count per sequence length
table <- as.data.frame(table(nchar(colnames(seqtab))))
colnames(table) <- c("LENGTH","COUNT")

ggplot(table,aes(x=LENGTH,y=COUNT)) + 
  geom_histogram(stat="identity") + 
  ggtitle("Sequence Lengths by SEQ Count") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=8)) +
  theme(axis.text.y=element_text(size=10))

# Plot total number of reads per sequence length
table2 <- tapply(colSums(seqtab), nchar(colnames(seqtab)), sum)
table2 <- data.frame(key=names(table2), value=table2)

colnames(table2) <- c("LENGTH","ABUNDANCE")

ggplot(table2,aes(x=LENGTH,y=ABUNDANCE)) + 
  geom_histogram(stat="identity") + 
  ggtitle("Sequence Lengths by SEQ Abundance") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=8)) +
  theme(axis.text.y=element_text(size=10))

##############################################################################
##### Remove columns with lengths outside the set sequence length window #####
##############################################################################
seqlens <- nchar(colnames(seqtab))
seqtab.filt <- seqtab[, (seqlens <= 254) & (seqlens >= 252)]
# Sanity check. Inspect distribution of sequence lengths
print(table(nchar(getSequences(seqtab.filt))))

###########################
##### Remove chimeras #####
###########################
seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab.filt)

############################################
##### Track reads through the pipeline #####
############################################
# As a final check, look at the number of reads that made it through each step in the pipeline.
getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# Raname columns
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
# Set rownames (use samples2 if you have to remove low-count samples)
rownames(track) <- samples
# View
print(track)
# Save
write.csv(track, "dada2_track_seqs.csv") 

############################################################################
##### SAVE THE FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS ##### 
############################################################################
saveRDS(seqtab.nochim, file="~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/seqtab.nochim.rds")
# RELOAD THE SAVED INFO FROM HERE (if you have closed the project):
seqtab.nochim <- readRDS("~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/seqtab.nochim.rds")


##########################################
##########################################
##### II. ASSIGN TAXONOMY WITH DADA2 #####
##########################################
##########################################
# Downloaded the Silva version 138 database here https://zenodo.org/record/3986799#.X6X7wFNKhBw
taxa <- assignTaxonomy(seqtab.nochim, "~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
# Make species level assignments based on exact matching between ASVs and sequenced reference strains. 
# Recent analysis suggests that exact matching (or 100% identity) is the only appropriate way to assign 
# species to 16S gene fragments. Currently, species-assignment training fastas are available for the Silva 
# and RDP 16S databases. Download the silva_species_assignment_v138.fa.gz file.
taxa <- addSpecies(taxa, "~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/silva_species_assignment_v138.fa.gz")

##########################################
##### Export tax table and asv table #####
##########################################
# Inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
# Export the file
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
# Export the files
write.table(taxon,"~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/silva_otu_table.txt",sep="\t",col.names=NA)

########################################################
##### Export a separate tax table with NAs removed #####
########################################################
# Fix the NAs in the taxa table
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]
taxon$Species[is.na(taxon$Species)] <- taxon$Genus[is.na(taxon$Species)]
# Inspect the taxonomic assignments again:
taxon.print <- taxon # Removing sequence rownames for display only
rownames(taxon.print) <- NULL
head(taxon.print)
# Export the file with NAs removed
write.table(taxon,"~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/silva_taxa_table_noNAS.txt",sep="\t",col.names=NA)


############################################
############################################
##### III. CONSTRUCT PHYLOGENETIC TREE #####
############################################
############################################
# https://compbiocore.github.io/metagenomics-workshop/assets/DADA2_tutorial.html
# https://github.com/benjjneb/dada2/issues/88

###############################################
##### Extract sequences from DADA2 output #####
###############################################
sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences

#######################################################
##### Run Sequence Alignment (MSA) using DECIPHER #####
#######################################################
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

####################################################################
##### Change sequence alignment output into a phyDat structure ##### 
####################################################################
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

##################################
##### Create distance matrix #####
##################################
dm <- dist.ml(phang.align)

####################################
##### Perform Neighbor joining #####
####################################
treeNJ <- NJ(dm) # Note, tip order != sequence order

#######################################
##### Internal maximum likelihood #####
#######################################
fit = pml(treeNJ, data=phang.align)
# negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))


############################################################################
##### SAVE THE FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS #####
############################################################################
saveRDS(fitGTR, file="~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/fitGTR.rds")
# RELOAD THE SAVED INFO FROM HERE (if you have closed the project):
fitGTR <- readRDS("~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/fitGTR.rds")


#####################################
#####################################
#### IV. TEST LOAD INTO PHYLOSEQ ####
#####################################
#####################################
# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/silva_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/silva_taxa_table_noNAS.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("~/Box/Kemp_MSEW_2020/new_analysis/Rcode/_bin/mapping_file.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon),phy_tree(fitGTR$tree))
