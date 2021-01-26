########################
#### LOAD Libraries ####
########################
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(compositions)
library(zCompositions)
library(CoDaSeq)


########################
########################
#### ANCOM FUNCTION ####
########################
########################

#https://github.com/sidhujyatha/ANCOM

ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

############################
############################
#### LOAD INTO PHYLOSEQ ####
############################
############################

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("_bin/silva_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("_bin/silva_taxa_table_noNAS.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("_bin/mapping_file.txt",sep="\t",header=T,row.names=1)
fitGTR <- readRDS("_bin/fitGTR.rds")
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon),phy_tree(fitGTR$tree))
ps

# add a DNAStringSet object (to be accessed by the refseq function) to keep track of the ASV sequences
sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)
ps

# rename our ASVs to something more convenient in downstream analysis, while automatically retaining the corresponding unique sequence identifier
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# remove chloroplasts, mitochondria, and Eukaryota
ps <- subset_taxa(ps, Kingdom !="NA")
ps <- subset_taxa(ps, Family !="Mitochondria")
ps <- subset_taxa(ps, Order !="Chloroplast")

# filter very low abundance ASVs
ps2<-filter_taxa(ps, function(x) mean(x) >5, TRUE)

# add read count
reads_sample <- microbiome::readcount(ps2)
sample_data(ps2)$reads_sample <- reads_sample
sample_data(ps2)

# remove samples with less than 3,000 reads
ps2 =subset_samples(ps2, reads_sample >3000)
min <- min(sample_sums(ps2))
sample_data(ps2) 

###################################
###################################
#### ANCOM at the family level ####
###################################
###################################

###########################
#### Family level PD28 ####
###########################
sub_ps = subset_samples(ps2, PD=="PD28")
fam = tax_glom(sub_ps , taxrank = "Family")
get_taxa_unique(sub_ps, "Family")
sub_ps <- subset_taxa(sub_ps, Family!="NA")
get_taxa_unique(sub_ps, "Family")

# Prepare count table
otu = as(otu_table(fam), "matrix")
otu_data = as.data.frame((otu))
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new dataframe for taxon
ID.taxa <- as.data.frame(tax_table(fam))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(fam), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=TRUE,
                           repeated=F,
                           main.var="TREATMENT",
                           adj.formula="MOM",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.7==TRUE),] 
res = res %>% dplyr::rename("taxa_id" = otu.names) 
res <- left_join(res,  ID.taxa, by="taxa_id")
pd28_ancom_family <- res
print(res2$Taxa)
write.table(pd28_ancom_family,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/PD28_ANCOM_family.txt",sep="\t",col.names=NA)
write.table(res,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/res_pd28_fam.txt",sep="\t",col.names=NA)

# Explore ASVs within significant Families
ps_pepto <- subset_taxa(ps2, Family=="Peptostreptococcaceae")
print(refseq(ps_pepto))

###########################
#### Family level PD10 ####
###########################
sub_ps = subset_samples(ps2, PD=="PD10")
fam = tax_glom(sub_ps , taxrank = "Family")
get_taxa_unique(sub_ps, "Family")
sub_ps <- subset_taxa(sub_ps, Family!="NA")
get_taxa_unique(sub_ps, "Family")

# Prepare count table
otu = as(otu_table(fam), "matrix")
otu_data = as.data.frame((otu))
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new dataframe for taxon
ID.taxa <- as.data.frame(tax_table(fam))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(fam), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=TRUE,
                           repeated=F,
                           main.var="TREATMENT",
                           adj.formula="MOM",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.7==TRUE),] 
res = res %>% dplyr::rename("taxa_id" = otu.names) 
res <- left_join(res,  ID.taxa, by="taxa_id")
pd10_ancom_family <- res
print(res2$Taxa)
write.table(pd10_ancom_family,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/PD10_ANCOM_family.txt",sep="\t",col.names=NA)
write.table(res,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/res_pd10_fam.txt",sep="\t",col.names=NA)

# Explore ASVs within significant Families
ps_entero <- subset_taxa(ps2, Family=="Enterococcaceae")
print(refseq(ps_entero))


##########################
#### Family level PD2 ####
##########################
sub_ps = subset_samples(ps2, PD=="PD2")
fam = tax_glom(sub_ps , taxrank = "Family")
get_taxa_unique(sub_ps, "Family")
sub_ps <- subset_taxa(sub_ps, Family!="NA")
get_taxa_unique(sub_ps, "Family")

# Prepare count table
otu = as(otu_table(fam), "matrix")
otu_data = as.data.frame((otu))
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new dataframe for taxon
ID.taxa <- as.data.frame(tax_table(fam))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(fam), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=TRUE,
                           repeated=F,
                           main.var="TREATMENT",
                           adj.formula="MOM",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.7==TRUE),] 
res2 = res2 %>% dplyr::rename("taxa_id" = otu.names) 
res2 <- left_join(res2,  ID.taxa, by="taxa_id")
print(res2$Taxa)
character(0)

###########################
#### Family level Dams ####
###########################

## Family level Dams - PAIRED TEST
sub_ps = subset_samples(ps2, PD=="Dam" & MOM!="dam_14")
fam = tax_glom(sub_ps , taxrank = "Family")
get_taxa_unique(sub_ps, "Family")
sub_ps <- subset_taxa(sub_ps, Family !="NA")
get_taxa_unique(sub_ps, "Family")

# Prepare count table
otu = as(otu_table(fam), "matrix")
otu_data = as.data.frame(otu)
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new dataframe for taxon
ID.taxa <- as.data.frame(tax_table(fam))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(fam), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=F,
                           repeated=T,
                           main.var="TREATMENT",
                           adj.formula=NULL,
                           repeat.var="MOM",
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.7==TRUE),] 
res2 = res2 %>% dplyr::rename("taxa_id" = otu.names) 
res2 <- left_join(res2,  ID.taxa, by="taxa_id")
print(res2$Taxa)
# character(0)


###################################
###################################
#### ANCOM at the family level ####
###################################
###################################

##########################
#### Genus level PD28 ####
##########################
sub_ps = subset_samples(ps2, PD=="PD28")
gen = tax_glom(sub_ps , taxrank = "Genus")
get_taxa_unique(sub_ps, "Genus")
sub_ps <- subset_taxa(sub_ps, Genus!="NA")
get_taxa_unique(sub_ps, "Genus")

# Prepare count table
otu = as(otu_table(gen), "matrix")
otu_data = as.data.frame((otu))
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new dataframe for taxon
ID.taxa <- as.data.frame(tax_table(gen))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(gen), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=TRUE,
                           repeated=F,
                           main.var="TREATMENT",
                           adj.formula="MOM",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.7==TRUE),] 
res = res %>% dplyr::rename("taxa_id" = otu.names) 
res <- left_join(res,  ID.taxa, by="taxa_id")
pd28_ancom_genus <- res
print(res2$Taxa)
write.table(pd28_ancom_genus,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/PD28_ANCOM_genus.txt",sep="\t",col.names=NA)
write.table(res,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/res_pd28_g.txt",sep="\t",col.names=NA)

# Explore ASVs of significant genera
ps_UCG005 <- subset_taxa(ps2, Genus=="UCG-005")
print(refseq(ps_UCG005))
ps_Tyz <- subset_taxa(ps2, Genus=="Tyzzerella")
refseq(ps_Tyz)
ps_Rumino <- subset_taxa(ps2, Genus=="Ruminococcus")
refseq(ps_Rumino)

##########################
#### Genus level PD10 ####
##########################
sub_ps = subset_samples(ps2, PD=="PD10")
gen = tax_glom(sub_ps , taxrank = "Genus")
get_taxa_unique(sub_ps, "Genus")
sub_ps <- subset_taxa(sub_ps, Genus !="NA")
get_taxa_unique(sub_ps, "Genus")

# Prepare count table
otu = as(otu_table(gen), "matrix")
otu_data = as.data.frame(otu)
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new dataframe for taxon
ID.taxa <- as.data.frame(tax_table(gen))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:7, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(gen), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=TRUE,
                           repeated=F,
                           main.var="TREATMENT",
                           adj.formula="MOM",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.7==TRUE),] 
res = res %>% dplyr::rename("taxa_id" = otu.names) 
res <- left_join(res,  ID.taxa, by="taxa_id")
pd10_ancom_genus <- res
print(res2$Taxa)
write.table(pd10_ancom_genus,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/PD10_ANCOM_genus.txt",sep="\t",col.names=NA)
write.table(res,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/res_pd10_g.txt",sep="\t",col.names=NA)

# Explore ASVs of significant genera
ps_Rumino <- subset_taxa(ps2, Genus=="Enterococcus")
refseq(ps_Rumino)


#########################
#### Genus level PD2 ####
#########################
sub_ps = subset_samples(ps2, PD=="PD2")
gen = tax_glom(sub_ps , taxrank = "Genus")
get_taxa_unique(sub_ps, "Genus")
sub_ps <- subset_taxa(sub_ps, Genus !="NA")
get_taxa_unique(sub_ps, "Genus")

# Prepare count table
otu = as(otu_table(gen), "matrix")
otu_data = as.data.frame(otu)
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new dataframe for taxon
ID.taxa <- as.data.frame(tax_table(gen))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(gen), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=TRUE,
                           repeated=F,
                           main.var="TREATMENT",
                           adj.formula="MOM",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.7==TRUE),] 
res2 = res2 %>% dplyr::rename("taxa_id" = otu.names) 
res2 <- left_join(res2,  ID.taxa, by="taxa_id")
print(res2$Taxa)
# character(0)

##########################
#### Genus level Dams ####
##########################

## Genus level Dams - PAIRED TEST
sub_ps = subset_samples(ps2, PD=="Dam" & MOM!="dam_14")
gen = tax_glom(sub_ps , taxrank = "Genus")
get_taxa_unique(sub_ps, "Genus")
sub_ps <- subset_taxa(sub_ps, Genus !="NA")
get_taxa_unique(sub_ps, "Genus")

# Prepare count table
otu = as(otu_table(gen), "matrix")
otu_data = as.data.frame(otu)
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new dataframe for taxon
ID.taxa <- as.data.frame(tax_table(gen))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(gen), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=F,
                           repeated=T,
                           main.var="TREATMENT",
                           adj.formula=NULL,
                           repeat.var="MOM",
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.7==TRUE),] 
res2 = res2 %>% dplyr::rename("taxa_id" = otu.names) 
res2 <- left_join(res2,  ID.taxa, by="taxa_id")
print(res2$Taxa)
# character(0)


################################
################################
#### ANCOM at the ASV level ####
################################
################################

########################
#### ASV level PD28 ####
########################
sub_ps = subset_samples(ps2, PD=="PD28")
sub_ps = metagMisc::phyloseq_filter_prevalence(sub_ps, prev.trh = 0.30)

# Prepare count table
otu = as(otu_table(sub_ps), "matrix")
otu_data = as.data.frame((otu))
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new datafram for taxon
ID.taxa <- as.data.frame(tax_table(sub_ps))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(sub_ps), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=TRUE,
                           repeated=F,
                           main.var="TREATMENT",
                           adj.formula="MOM",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
write.table(res,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/res_pd28_asv.txt",sep="\t",col.names=NA)
res2 <- res
res2 <- res[which(res$detected_0.7==TRUE),] 
res2 = res2 %>% dplyr::rename("taxa_id" = otu.names) 
res2 <- left_join(res2,  ID.taxa, by="taxa_id")
pd28_ancom_asv <- res2
pd28_asv_list <- res2 %>% pull(taxa_id)
print(res2$Taxa)
write.table(pd28_ancom_asv,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/PD28_ANCOM_asv.txt",sep="\t",col.names=NA)

########################
#### ASV level PD10 ####
########################
sub_ps = subset_samples(ps2, PD=="PD10")
sub_ps = metagMisc::phyloseq_filter_prevalence(sub_ps, prev.trh = 0.30)

# Prepare count table
otu = as(otu_table(sub_ps), "matrix")
otu_data = as.data.frame((otu))
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new datafram for taxon
ID.taxa <- as.data.frame(tax_table(sub_ps))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(sub_ps), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=TRUE,
                           repeated=F,
                           main.var="TREATMENT",
                           adj.formula="MOM",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.6==TRUE),] 
res2 = res2 %>% dplyr::rename("taxa_id" = otu.names) 
res2 <- left_join(res2,  ID.taxa, by="taxa_id")
pd10_ancom_asv <- res2
pd10_asv_list <- pd10_ancom_asv %>% pull(taxa_id)
print(res2$Taxa)
write.table(pd10_ancom_asv,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/PD10_ANCOM_asv.txt",sep="\t",col.names=NA)

#######################
#### ASV level PD2 ####
#######################
sub_ps = subset_samples(ps2, PD=="PD2")
sub_ps = metagMisc::phyloseq_filter_prevalence(sub_ps, prev.trh = 0.30)

# Prepare count table
otu = as(otu_table(sub_ps), "matrix")
otu_data = as.data.frame((otu))
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new datafram for taxon
ID.taxa <- as.data.frame(tax_table(sub_ps))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(sub_ps), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=TRUE,
                           repeated=F,
                           main.var="TREATMENT",
                           adj.formula="MOM",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.6==TRUE),] 
res2 = res2 %>% dplyr::rename("taxa_id" = otu.names) 
res2 <- left_join(res2,  ID.taxa, by="taxa_id")
print(res2$Taxa)
# character(0)

########################
#### ASV level Dams ####
########################

## ASV level Dam - PAIRED TEST
sub_ps = subset_samples(ps2, PD=="Dam" & MOM!="dam_14") #remove dam that had only one litter for paired test
sub_ps = metagMisc::phyloseq_filter_prevalence(sub_ps, prev.trh = 0.30)

# Prepare count table
otu = as(otu_table(sub_ps), "matrix")
otu_data = as.data.frame((otu))
# make row names a column with the heading "Sample.ID)
otu_data <- otu_data %>% tibble::rownames_to_column("Sample.ID")

# Make new datafram for taxon
ID.taxa <- as.data.frame(tax_table(sub_ps))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 3:8, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)

# Prepare metadata
metadata = as(sample_data(sub_ps), "matrix")
meta_data= as.data.frame(metadata)
meta_data = meta_data %>% dplyr::rename("Sample.ID" = SampleID) %>% 
  dplyr::select(1:6)
rownames(meta_data) <- c()

# Run ANCOM
comparison_test=ANCOM.main(OTUdat=otu_data,
                           Vardat=meta_data,
                           adjusted=F,
                           repeated=T,
                           main.var="TREATMENT",
                           adj.formula=NULL,
                           repeat.var="MOM",
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

res <- comparison_test$W.taxa #taxa that significantly vary across factor level of interest
res2 <- res[which(res$detected_0.7==TRUE),] 
res2 = res2 %>% dplyr::rename("taxa_id" = otu.names) 
res2 <- left_join(res2,  ID.taxa, by="taxa_id")
print(res2$Taxa)


###############################
#### Heatmap for PD28 ASVs ####
###############################

#CLT transformation
ps_pd28 = subset_samples(ps2, PD=="PD28")
pd28_otu = as(otu_table(ps_pd28), "matrix")
otu = as.data.frame(pd28_otu)
metadata = as.matrix(sample_data(ps_pd28))
meta_data = as.data.frame(metadata)
# count number of zeros in the table
sum(otu == 0)
#[1] 3402
# Replace 0 values with an estimate (because normalization is taking log, can't have 0) using the Bayesian-multiplicative replacement function cmultRepl in Zcompositions. Also transposing here, because we need samples as rows.
d.czm <- zCompositions::cmultRepl(t(otu), method="CZM", label=0)
#No. corrected values:  2701 
# Perform the centred log-ratio, or CLR (Aitchison) transformation using CoDaSeq.  
# The output will have the ASVs as columns and samples as ROWS
d.clr <- CoDaSeq::codaSeq.clr(d.czm)
dim(d.clr)
#[1] 370  25
# taxa as row samples as columns

# Convert transformed object to a data frame
df.clr <- as.data.frame(d.clr)
# Copy rownames to a column
df.clr$taxa_id<- rownames(df.clr)
head(df.clr)

# Keep only the significantly differentiated ASVs at PD28
sig_ASVs_pd28 <- df.clr[df.clr$taxa_id %in% pd28_asv_list,]
mat2 <- sig_ASVs_pd28 

# get annotations
ID.taxa <- as.data.frame(tax_table(ps2))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 6:7, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)
ID.taxa_pd28 <- ID.taxa[ID.taxa$taxa_id %in% pd28_asv_list,]
ID.taxa_pd28$taxa_id <- ordered(ID.taxa_pd28$taxa_id , levels=pd28_asv_list)
mat2$taxa_id <- ordered(mat2$taxa_id , levels=pd28_asv_list)
# check to make sure order matches
ID.taxa_pd28$taxa_id == mat2$taxa_id
# should be TRUE for all

# row annotations
annotation_row <- ID.taxa[ID.taxa$taxa_id %in% pd28_asv_list,]
row.names(annotation_row) <- NULL
annotation_row <- annotation_row %>% column_to_rownames("taxa_id")
# set row names in mat2
mat2$taxa_id <- NULL
annotation_row <- subset(annotation_row, select = Taxa)
unique(annotation_row$Taxa)

# col annotations
annotation_col <- meta_data %>% dplyr::select(TREATMENT)

# https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
# assign colors
cols_rows = c( "Muribaculaceae, Muribaculaceae" = "#696969",                        
               "Lachnospiraceae, Lachnospiraceae UCG-001" = "#2e8b57",             
               "Lachnospiraceae, [Eubacterium] xylanophilum group" =  "#7fff00" , 
               "Lachnospiraceae, Lachnospiraceae NK4A136 group"   =    "#00bfff" ,
               "Lachnospiraceae, Lachnoclostridium"     =    "#0000cd" ,
               "Lachnospiraceae, Lachnospiraceae"  =    "#00fa9a",
               "Lachnospiraceae, Marvinbryantia"  =  "#4169e1",                    
               "Lachnospiraceae, Tyzzerella"   =  "#00ffff" , 
               "Ruminococcaceae, Incertae Sedis"   =   "#ffa500",  
               "Ruminococcaceae, Ruminococcaceae" =     "#ffff00",
               "Oscillospiraceae, UCG-005"  =    "#ff00ff",
               "Oscillospiraceae, Oscillibacter"  =   "#800080", 
               "Oscillospiraceae, Intestinimonas"  =   "#dda0dd",
               "Clostridiaceae, Clostridium sensu stricto 1" =   "#ff1493" ,      
               "Clostridia vadinBB60 group, Clostridia vadinBB60 group" = "#8b0000",
               "Clostridia UCG-014, Clostridia UCG-014"  =  "#e9967a"              
)                
my_colour = list(
  Taxa = cols_rows,
  TREATMENT = c("Control" = "black", "MSEW" = "red")
)

library(pheatmap)

hm2 <- pheatmap(mat2,annotation_row = annotation_row,annotation_col = annotation_col,
                clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", border_color = "grey60",
                clustering_method = 'complete', fontsize =10, cluster_cols = T, cluster_rows = T,annotation_colors = my_colour,
                cellwidth = 20, cutree_rows = 2, cutree_cols = 2, cellheight = 10)
hm2


###############################
#### Heatmap for PD10 ASVs ####
###############################

#CLT transformation
ps_pd10 = subset_samples(ps2, PD=="PD10")
pd10_otu = as(otu_table(ps_pd10), "matrix")
otu = as.data.frame(pd10_otu)
metadata = as.matrix(sample_data(ps_pd10))
meta_data = as.data.frame(metadata)
# count number of zeros in the table
sum(otu == 0)
#[1] 3402
# Replace 0 values with an estimate (because normalization is taking log, can't have 0) using the Bayesian-multiplicative replacement function cmultRepl in Zcompositions. Also transposing here, because we need samples as rows.
d.czm <- zCompositions::cmultRepl(t(otu), method="CZM", label=0)
#No. corrected values:  2701 
# Perform the centred log-ratio, or CLR (Aitchison) transformation using CoDaSeq.  
# The output will have the ASVs as columns and samples as ROWS
d.clr <- CoDaSeq::codaSeq.clr(d.czm)
dim(d.clr)
#[1] 370  25
# taxa as row samples as columns

# Convert transformed object to a data frame
df.clr <- as.data.frame(d.clr)
# Copy rownames to a column
df.clr$taxa_id<- rownames(df.clr)
head(df.clr)


# Keep only the significantly differentiated ASVs at PD28
sig_ASVs_pd10 <- df.clr[df.clr$taxa_id %in% pd10_asv_list,]
mat2 <- sig_ASVs_pd10 

# get annotations
ID.taxa <- as.data.frame(tax_table(ps2))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 6:7, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)
ID.taxa_pd10 <- ID.taxa[ID.taxa$taxa_id %in% pd10_asv_list,]
ID.taxa_pd10$taxa_id <- ordered(ID.taxa_pd10$taxa_id , levels=pd10_asv_list)
mat2$taxa_id <- ordered(mat2$taxa_id , levels=pd10_asv_list)
# check to make sure order matches
ID.taxa_pd10$taxa_id == mat2$taxa_id
# should be TRUE for all

# row annotation
annotation_row <- ID.taxa[ID.taxa$taxa_id %in% pd10_asv_list,]
row.names(annotation_row) <- NULL
annotation_row <- annotation_row %>% column_to_rownames("taxa_id")
# set row names in mat2
mat2$taxa_id <- NULL
annotation_row <- subset(annotation_row, select = Taxa)
unique(annotation_row$Taxa)

# col annotations
annotation_col <- meta_data %>% dplyr::select(TREATMENT)
#annotation_col <- annotation_col %>% column_to_rownames('SampleID')


# https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
# assign colors
cols_rows = c( "Enterococcaceae, Enterococcus" = "#eee8aa"                      
)                
my_colour = list(
  Taxa = cols_rows,
  TREATMENT = c("Control" = "black", "MSEW" = "red")
)

hm3 <- pheatmap(mat2,annotation_row = annotation_row,annotation_col = annotation_col,
                clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", border_color = "grey60",
                clustering_method = 'complete', fontsize =10, cluster_cols = T, cluster_rows = F,annotation_colors = my_colour,
                cellwidth = 20, cutree_rows = 1, cutree_cols = 3, cellheight = 10)
hm3


###############################
#### Heatmap for Dam ASVs ####
###############################

# make a chart of "PD10" and "PD28" OTUs in Dam samples
sigASVs <- c(pd10_asv_list,pd28_asv_list)

#CLT transformation
ps_dam = subset_samples(ps2, PD=="Dam")
dam_otu = as(otu_table(ps_dam), "matrix")
otu = as.data.frame(dam_otu)
metadata = as.matrix(sample_data(ps_dam))
meta_data = as.data.frame(metadata)
# count number of zeros in the table
sum(otu == 0)
#[1] 3402
# Replace 0 values with an estimate (because normalization is taking log, can't have 0) using the Bayesian-multiplicative replacement function cmultRepl in Zcompositions. Also transposing here, because we need samples as rows.
d.czm <- zCompositions::cmultRepl(t(otu), method="CZM", label=0)
#No. corrected values:  2701 
# Perform the centred log-ratio, or CLR (Aitchison) transformation using CoDaSeq.  
# The output will have the ASVs as columns and samples as ROWS
d.clr <- CoDaSeq::codaSeq.clr(d.czm)
dim(d.clr)
#[1] 370  25
# taxa as row samples as columns

# Convert transformed object to a data frame
df.clr <- as.data.frame(d.clr)
# Copy rownames to a column
df.clr$taxa_id<- rownames(df.clr)
head(df.clr)

# Keep only the significantly differentiated ASVs at PD28 and PD10
sig_ASVs <- df.clr[df.clr$taxa_id %in% sigASVs,]
mat2 <- sig_ASVs

# get annotations
ID.taxa <- as.data.frame(tax_table(ps2))
ID.taxa<- rownames_to_column(ID.taxa, var = "taxa_id") 
ID.taxa <- ID.taxa %>% unite("Taxa", 6:7, sep = ", ", remove = FALSE) %>% dplyr::select(taxa_id,Taxa)
ID.taxa_dam <- ID.taxa[ID.taxa$taxa_id %in% sigASVs,]
ID.taxa_dam$taxa_id <- ordered(ID.taxa_dam$taxa_id , levels=sigASVs)
mat2$taxa_id <- ordered(mat2$taxa_id , levels=sigASVs)
# check to make sure order matches
ID.taxa_dam$taxa_id == mat2$taxa_id
# should be TRUE for all

# row annotation
annotation_row <- ID.taxa[ID.taxa$taxa_id %in% sigASVs,]
row.names(annotation_row) <- NULL
annotation_row <- annotation_row %>% column_to_rownames("taxa_id")
# set row names in mat2
mat2$taxa_id <- NULL

# set order
col_order <- c("MOM-011-2-2018","MOM-012-2-2018","MOM-013-2-2018","Mom-015-5-2-18","Mom-016-5-2-18",
               "Mom-011-5-23-18","Mom-012-5-23-18","Mom-013-5-24-18","MOM-015-2-2018","MOM-016-2-2018")
mat2 <- mat2[, col_order]


annotation_row <- subset(annotation_row, select = Taxa)
unique(annotation_row$Taxa)
write.csv(annotation_row,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/dam_annotation_row.csv",row.names=TRUE)
annotation_row2 <- read.csv("~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/dam_annotation_row.csv")
annotation_row2 <- annotation_row2 %>% column_to_rownames("X")

# col annotations
annotation_col <- meta_data %>% dplyr::select(TREATMENT)


# https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
# assign colors
cols_rows = c( "Muribaculaceae, Muribaculaceae" = "#696969",                        
               "Lachnospiraceae, Lachnospiraceae UCG-001" = "#2e8b57",             
               "Lachnospiraceae, [Eubacterium] xylanophilum group" =  "#7fff00" , 
               "Lachnospiraceae, Lachnospiraceae NK4A136 group"   =    "#00bfff" ,
               "Lachnospiraceae, Lachnoclostridium"     =    "#0000cd" ,
               "Lachnospiraceae, Lachnospiraceae"  =    "#00fa9a",
               "Lachnospiraceae, Marvinbryantia"  =  "#4169e1",                    
               "Lachnospiraceae, Tyzzerella"   =  "#00ffff" , 
               "Ruminococcaceae, Incertae Sedis"   =   "#ffa500",  
               "Ruminococcaceae, Ruminococcaceae" =     "#ffff00",
               "Oscillospiraceae, UCG-005"  =    "#ff00ff",
               "Oscillospiraceae, Oscillibacter"  =   "#800080", 
               "Oscillospiraceae, Intestinimonas"  =   "#dda0dd",
               "Clostridiaceae, Clostridium sensu stricto 1" =   "#ff1493" ,      
               "Clostridia vadinBB60 group, Clostridia vadinBB60 group" = "#8b0000",
               "Clostridia UCG-014, Clostridia UCG-014"  =  "#e9967a",
               "Enterococcaceae, Enterococcus" = "#eee8aa"
) 
my_colour = list(
  Trend.in.MSEW.Pups = c("increased abundance" = "red", "decreased abundance" = "black"),
  Taxa = cols_rows,
  TREATMENT = c("Control" = "grey", "MSEW" = "#A32CC4")
)

# define col labels
labels_col <- c("Dam 011","Dam 012","Dam 013",
                "Dam 015", "Dam 016",
                "Dam 011", "Dam 012","Dam 013",
                "Dam 015", "Dam 016")

hm3 <- pheatmap(mat2,annotation_row = annotation_row2, annotation_col = annotation_col,
                clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", border_color = "grey60",
                clustering_method = 'complete', fontsize =10, cluster_cols = F, cluster_rows = T,annotation_colors = my_colour,
                cellwidth = 20, labels_col = labels_col, cutree_rows = 1, cutree_cols = 1, cellheight = 10)
hm3


###############################
#### Volcano plot for ASVs ####
###############################

#CLT transformation
ps_pd28 = subset_samples(ps2, PD=="PD28")

pd28_asv = as(otu_table(ps_pd28), "matrix")
otu = as.data.frame(pd28_asv)
metadata = as.matrix(sample_data(ps_pd28))
meta_data = as.data.frame(metadata)
# count number of zeros in the table
sum(otu == 0)
#[1] 3402
# Replace 0 values with an estimate (because normalization is taking log, can't have 0) using the Bayesian-multiplicative replacement function cmultRepl in Zcompositions. Also transposing here, because we need samples as rows.
d.czm <- zCompositions::cmultRepl(t(otu), method="CZM", label=0)
#No. corrected values:  2701 
# Perform the centred log-ratio, or CLR (Aitchison) transformation using CoDaSeq.  
# The output will have the ASVs as columns and samples as ROWS
d.clr <- CoDaSeq::codaSeq.clr(d.czm)
dim(d.clr)
#[1] 370  25
# taxa as row samples as columns

# Convert transformed object to a data frame
df.clr <- as.data.frame(d.clr)
# Copy rownames to a column
df.clr$taxa_id<- rownames(df.clr)
head(df.clr)
write_csv(df.clr,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/ANCOM_pd28_asv.csv")
vplot_pd28_asv<-read_csv("~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/res_pd28_asv_volcano_plot_upload.csv")

theme_set(theme_bw())
cols = c("FALSE" = "grey", "TRUE" = "red")
pd28_volcano_asv <- ggplot(vplot_pd28_asv) +
  geom_point(size=3,aes(x=clr, y=W_stat, colour=detected_0.7)) +
  ggtitle("") +
  xlab("Centered Log-Ratio (CLR) Mean Difference") + 
  ylab("W Statistic") +
  scale_colour_manual(values=cols)+
  theme(legend.position = "none") +
  scale_x_continuous(limits=c(-5,5), breaks = seq(-5, 5, by = 1))
pd28_volcano_asv


######################################
#### Volcano plot for PD28 Family ####
######################################

#CLT transformation
ps_pd28 = subset_samples(ps2, PD=="PD28")
fam = tax_glom(ps_pd28 , taxrank = "Family")
get_taxa_unique(fam, "Family")


pd28_fam = as(otu_table(fam), "matrix")
otu = as.data.frame(pd28_fam)
metadata = as.matrix(sample_data(ps_pd28))
meta_data = as.data.frame(metadata)
# count number of zeros in the table
sum(otu == 0)
#[1] 3402
# Replace 0 values with an estimate (because normalization is taking log, can't have 0) using the Bayesian-multiplicative replacement function cmultRepl in Zcompositions. Also transposing here, because we need samples as rows.
d.czm <- zCompositions::cmultRepl(t(otu), method="CZM", label=0)
#No. corrected values:  2701 
# Perform the centred log-ratio, or CLR (Aitchison) transformation using CoDaSeq.  
# The output will have the ASVs as columns and samples as ROWS
d.clr <- CoDaSeq::codaSeq.clr(d.czm)
dim(d.clr)
#[1] 370  25
# taxa as row samples as columns

# Convert transformed object to a data frame
df.clr <- as.data.frame(d.clr)
# Copy rownames to a column
df.clr$taxa_id<- rownames(df.clr)
head(df.clr)
write_csv(df.clr,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/ANCOM_pd28_fam.csv")
vplot_pd28_fam<-read_csv("~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/res_pd28_fam_volcano_plot_upload.csv")

theme_set(theme_bw())
cols = c("FALSE" = "grey", "TRUE" = "red")
pd28_volcano_fam <- ggplot(vplot_pd28_fam) +
  geom_point(size=3,aes(x=crl, y=W_stat, colour=detected_0.7)) +
  ggtitle("") +
  xlab("Centered Log-Ratio (CLR) Mean Difference") + 
  ylab("W Statistic") +
  scale_colour_manual(values=cols)+
  theme(legend.position = "none") +
  scale_x_continuous(limits=c(-2.25,3), breaks = seq(-2, 3, by = 1))
pd28_volcano_fam
# export 650 x 450


#####################################
#### Volcano plot for PD28 Genus ####
#####################################

ps_pd28 = subset_samples(ps2, PD=="PD28")
gen = tax_glom(ps_pd28 , taxrank = "Genus")
get_taxa_unique(gen, "Genus")

pd28_gen = as(otu_table(gen), "matrix")
otu = as.data.frame(pd28_gen)
metadata = as.matrix(sample_data(ps_pd28))
meta_data = as.data.frame(metadata)
# count number of zeros in the table
sum(otu == 0)
#[1] 3402
# Replace 0 values with an estimate (because normalization is taking log, can't have 0) using the Bayesian-multiplicative replacement function cmultRepl in Zcompositions. Also transposing here, because we need samples as rows.
d.czm <- zCompositions::cmultRepl(t(otu), method="CZM", label=0)
#No. corrected values:  2701 
# Perform the centred log-ratio, or CLR (Aitchison) transformation using CoDaSeq.  
# The output will have the ASVs as columns and samples as ROWS
d.clr <- CoDaSeq::codaSeq.clr(d.czm)
dim(d.clr)
#[1] 370  25
# taxa as row samples as columns

# Convert transformed object to a data frame
df.clr <- as.data.frame(d.clr)
# Copy rownames to a column
df.clr$taxa_id<- rownames(df.clr)
head(df.clr)
write_csv(df.clr,"~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/ANCOM_pd28_gen.csv")
vplot_pd28<-read_csv("~/Box/Kemp_MSEW_2020/new_analysis/R_analysis/res_pd28_genus_volcano_plot_input.csv")

theme_set(theme_bw())
cols = c("FALSE" = "grey", "TRUE" = "black")
pd28_volcano <- ggplot(vplot_pd28) +
  geom_point(size=3,aes(x=clr, y=W_stat, colour=detected_0.7)) +
  ggtitle("") +
  xlab("Centered Log-Ratio (CLR) Mean Difference") + 
  ylab("W Statistic") +
  scale_colour_manual(values=cols)+
  theme(legend.position = "none") +
  scale_x_continuous(limits=c(-3.5,2), breaks = seq(-3, 2, by = 1))
pd28_volcano

