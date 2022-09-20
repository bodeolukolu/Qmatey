#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
libdir <- args[1]
.libPaths( c( .libPaths(), libdir) )
library(dplyr, quietly = T)
library(ggplot2, quietly = T)
library(reshape2, quietly = T)

#########################################################
# Visualization (sensitivity, false positive, and false negative) before cross-taxon validation
#########################################################
phylum <- read.delim(file="./phylum_level/phylum_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
class <- read.delim(file="./class_level/class_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
order <- read.delim(file="./order_level/order_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
family <- read.delim(file="./family_level/family_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
genus <- read.delim(file="./genus_level/genus_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
species <- read.delim(file="./species_level/soecies_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
strain <- read.delim(file="./strain_level_minUniq_1/strain_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
ctv <- as.data.frame(matrix(nrow=6,ncol=10))
colnames(ctv) <- c("ref","ref_len","subject","subject_len","overlap","ref_uniq","subjec_uniq","false_positive","false_negative","sensitivity")

taxalevels <- c("strain","species","genus","family","order","class","phylum")
for (i in c(1:6)){
  ref <- taxalevels[i+1]
  ref_len <- cat(taxalevels[i+1],taxalevels[i+1],sep="$"); ref_len <- length(ref_len[!is.na(ref_len)])
  subject <- taxalevels[i]
  subject_len <- cat(taxalevels[i],taxalevels[i],sep="$"); subject_len <- length(subject_len[!is.na(subject_len)])
  overlap <- intersect(cat(taxalevels[i+1],taxalevels[i+1],sep="$"), cat(taxalevels[i],taxalevels[i+1],sep="$")); overlap <- overlap(!is.na(overlap))
  ref_uniq <- !c(cat(taxalevels[i+1],taxalevels[i+1],sep="$")) %in% c(overlap); ref_uniq <- length(ref_uniq)
  subject_uniq <- !c(cat(taxalevels[i],taxalevels[i],sep="$")) %in% c(overlap); subject_uniq <- length(subject_uniq)
  overlap <- length(overlap)
  false_positive <- subject_uniq / subject_len
  false_negative <- ref_uniq / ref_len
  sensitivity <- overlap / (overlap + ref_uniq)
  ctv_hold <- c(ref,ref_len,subject,subject_len,overlap,ref_uniq,subject_uniq,sensitivity)
  ctv <- rbind(ctv,ctv_hold)
}
ctv$Cross_taxon_comparison <- paste(ctv$subject,"_vs._",ctv$ref,"_(ref)",sep="") 
ctv$Cross_taxon_comparison <- factor(a.character(ctv$Cross_taxon_comparison), levels=unique(ctv$Cross_taxon_comparison))
write.table(ctv,"cross_taxon_comparison_unvalidated.txt", sep="\t",row.names=FALSE, col.names=T, quote = F)
ctv <- subset(ctv, select=c(Cross_taxon_comparison,false_positive,false_negative,sensitivity))
ctv <- melt(ctv, id=Cross_taxon_comparison)

plot <- ggplot(ctv,aes(x=Cross_taxon_validation,y=value,group=variable,color=variable))+
  geom_line(size=1.5)+ylim(0,1)+
  geom_point()+xlab("Cross Taxon Comparison")+
  ylab("Value")+ theme_bw() + theme(legend.position="bottom") + 
  theme(text = element_text(size = 20),axis.text.x = element_text(angle=30,size=18.5,hjust = 1)) + labs(color='Metric')
ggsave(file="cross_taxon_comparison_unvalidated.tiff", plot=plot, width=12, height=6, units=("in"), dpi=300, compression = "lzw")


#########################################################
# Visualization (sensitivity, false positive, and false negative) after cross-taxon validation
#########################################################
phylum <- read.delim(file="./phylum_level_validated/phylum_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
class <- read.delim(file="./class_level_validated/class_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
order <- read.delim(file="./order_level_validated/order_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
family <- read.delim(file="./family_level_validated/family_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
genus <- read.delim(file="./genus_level_validated/genus_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
species <- read.delim(file="./species_level_validated/soecies_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
strain <- read.delim(file="./strain_level_minUniq_1_validated/strain_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
ctv <- as.data.frame(matrix(nrow=6,ncol=10))
colnames(ctv) <- c("ref","ref_len","subject","subject_len","overlap","ref_uniq","subjec_uniq","false_positive","false_negative","sensitivity")

taxalevels <- c("strain","species","genus","family","order","class","phylum")
for (i in c(1:6)){
  ref <- taxalevels[i+1]
  ref_len <- cat(taxalevels[i+1],taxalevels[i+1],sep="$"); ref_len <- length(ref_len[!is.na(ref_len)])
  subject <- taxalevels[i]
  subject_len <- cat(taxalevels[i],taxalevels[i],sep="$"); subject_len <- length(subject_len[!is.na(subject_len)])
  overlap <- intersect(cat(taxalevels[i+1],taxalevels[i+1],sep="$"), cat(taxalevels[i],taxalevels[i+1],sep="$")); overlap <- overlap(!is.na(overlap))
  ref_uniq <- !c(cat(taxalevels[i+1],taxalevels[i+1],sep="$")) %in% c(overlap); ref_uniq <- length(ref_uniq)
  subject_uniq <- !c(cat(taxalevels[i],taxalevels[i],sep="$")) %in% c(overlap); subject_uniq <- length(subject_uniq)
  overlap <- length(overlap)
  false_positive <- subject_uniq / subject_len
  false_negative <- ref_uniq / ref_len
  sensitivity <- overlap / (overlap + ref_uniq)
  ctv_hold <- c(ref,ref_len,subject,subject_len,overlap,ref_uniq,subject_uniq,sensitivity)
  ctv <- rbind(ctv,ctv_hold)
}
ctv$Cross_taxon_comparison <- paste(ctv$subject,"_vs._",ctv$ref,"_(ref)",sep="") 
ctv$Cross_taxon_comparison <- factor(a.character(ctv$Cross_taxon_comparison), levels=unique(ctv$Cross_taxon_comparison))
write.table(ctv,"cross_taxon_comparison_validated.txt", sep="\t",row.names=FALSE, col.names=T, quote = F)
ctv <- subset(ctv, select=c(Cross_taxon_comparison,false_positive,false_negative,sensitivity))
ctv <- melt(ctv, id=Cross_taxon_comparison)

plot <- ggplot(ctv,aes(x=Cross_taxon_validation,y=value,group=variable,color=variable))+
  geom_line(size=1.5)+ylim(0,1)+
  geom_point()+xlab("Cross Taxon Comparison")+
  ylab("Value")+ theme_bw() + theme(legend.position="bottom") + 
  theme(text = element_text(size = 20),axis.text.x = element_text(angle=30,size=18.5,hjust = 1)) + labs(color='Metric')
ggsave(file="cross_taxon_comparison_validated.tiff", plot=plot, width=12, height=6, units=("in"), dpi=300, compression = "lzw")

