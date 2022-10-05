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
species <- read.delim(file="./species_level/species_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
strain <- read.delim(file="./strain_level_minUniq_1/strain_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
phylum <- subset(phylum, domain != "Viruses")
class <- subset(class, domain != "Viruses")
order <- subset(order, domain != "Viruses")
family <- subset(family, domain != "Viruses")
genus <- subset(genus, domain != "Viruses")
species <- subset(species, domain != "Viruses")
strain <- subset(strain, domain != "Viruses")
ctv <- as.data.frame(matrix(nrow=0,ncol=10))

n <- max(length(phylum$phylum), length(class$class), length(order$order), length(family$family), length(genus$genus), length(species$species), length(strain$taxaname))
phylum <- phylum$phylum; length(phylum) <- n
class1 <- class$class; length(class1) <- n; class2 <- class$phylum; length(class2) <- n; class <- cbind(class1,class2)
order1 <- order$order; length(order1) <- n; order2 <- order$class; length(order2) <- n; order <- cbind(order1,order2)
family1 <- family$family; length(family1) <- n; family2 <- family$order; length(family2) <- n;  family <- cbind(family1,family2)
genus1 <- genus$genus; length(genus1) <- n; genus2 <- genus$family; length(genus2) <- n; genus <- cbind(genus1,genus2)
species1 <- species$species; length(species1) <- n; species2 <- species$genus; length(species2) <- n; species <- cbind(species1,species2)
strain <- strain$species; length(strain) <- n

taxalevels <- as.data.frame(cbind(strain,species,genus,family,order,class,phylum))
colnames(taxalevels) <- c("strain$species","species$species","species$genus","genus$genus","genus$family","family$family","family$order","order$order","order$class","class$class","class$phylum","phylum$phylum")
for (i in c(1,3,5,7,9,11)){
  ref <- gsub("\\$.*","",colnames(taxalevels)[i+1])
  ref_len <- taxalevels[,i+1]; ref_len <- length(ref_len[!is.na(ref_len)])
  subject <- gsub("\\$.*","",colnames(taxalevels)[i])
  subject_len <- taxalevels[,i]; subject_len <- length(subject_len[!is.na(subject_len)])
  overlap <- intersect(taxalevels[,i+1], taxalevels[,i]); overlap <- overlap[!is.na(overlap)]
  ref_uniq <- taxalevels[,i+1][!c(taxalevels[,i+1]) %in% c(overlap)]; ref_uniq <- length(ref_uniq[!is.na(ref_uniq)])
  subject_uniq <- taxalevels[,i][!c(taxalevels[,i]) %in% c(overlap)]; subject_uniq <- length(subject_uniq[!is.na(subject_uniq)])
  overlap <- length(overlap)
  false_positive <- subject_uniq / subject_len
  false_negative <- ref_uniq / ref_len
  sensitivity <- overlap / (overlap + ref_uniq)
  ctv_hold <- c(ref,ref_len,subject,subject_len,overlap,ref_uniq,subject_uniq,false_positive,false_negative,sensitivity)
  ctv <- rbind(ctv,ctv_hold)
}
colnames(ctv) <- c("ref","ref_len","subject","subject_len","overlap","ref_uniq","subjec_uniq","false_positive","false_negative","sensitivity")


ctv$Cross_taxon_comparison <- paste(ctv$subject,"_vs._",ctv$ref,"_(ref)",sep="") 
ctv$Cross_taxon_comparison <- factor(as.character(ctv$Cross_taxon_comparison), levels=unique(ctv$Cross_taxon_comparison))
write.table(ctv,"cross_taxon_comparison_unvalidated.txt", sep="\t",row.names=FALSE, col.names=T, quote = F)
ctv <- subset(ctv, select=c(Cross_taxon_comparison,false_positive,false_negative,sensitivity))
ctv <- melt(ctv, id="Cross_taxon_comparison")
ctv$value <- as.numeric(ctv$value)

plot <- ggplot(ctv,aes(x=Cross_taxon_comparison,y=value,group=variable,color=variable))+
  geom_line(size=1.5)+ylim(0,1)+
  geom_point()+xlab("Cross Taxon Comparison (unvalidated)")+
  ylab("Value")+ theme_bw() + theme(legend.position="bottom") + 
  theme(text = element_text(size = 20),axis.text.x = element_text(angle=30,size=18.5,hjust = 1)) + labs(color='Metric') +
  theme(plot.margin = margin(2,2,2,2, "cm"))
ggsave(file="cross_taxon_comparison_unvalidated.tiff", plot=plot, width=12, height=6, units=("in"), dpi=300, compression = "lzw")


#########################################################
# Visualization (sensitivity, false positive, and false negative) after cross-taxon validation
#########################################################
phylum <- read.delim(file="./phylum_level/phylum_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
class <- read.delim(file="./class_level_validated/class_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
order <- read.delim(file="./order_level_validated/order_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
family <- read.delim(file="./family_level_validated/family_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
genus <- read.delim(file="./genus_level_validated/genus_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
species <- read.delim(file="./species_level_validated/species_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
strain <- read.delim(file="./strain_level_minUniq_1_validated/strain_taxainfo_unique_sequences.txt", header=T, sep="\t", fill= T, quote="", check.names = T)
phylum <- subset(phylum, domain != "Viruses")
class <- subset(class, domain != "Viruses")
order <- subset(order, domain != "Viruses")
family <- subset(family, domain != "Viruses")
genus <- subset(genus, domain != "Viruses")
species <- subset(species, domain != "Viruses")
strain <- subset(strain, domain != "Viruses")
ctv <- as.data.frame(matrix(nrow=0,ncol=10))

n <- max(length(phylum$phylum), length(class$class), length(order$order), length(family$family), length(genus$genus), length(species$species), length(strain$taxaname))
phylum <- phylum$phylum; length(phylum) <- n
class1 <- class$class; length(class1) <- n; class2 <- class$phylum; length(class2) <- n; class <- cbind(class1,class2)
order1 <- order$order; length(order1) <- n; order2 <- order$class; length(order2) <- n; order <- cbind(order1,order2)
family1 <- family$family; length(family1) <- n; family2 <- family$order; length(family2) <- n;  family <- cbind(family1,family2)
genus1 <- genus$genus; length(genus1) <- n; genus2 <- genus$family; length(genus2) <- n; genus <- cbind(genus1,genus2)
species1 <- species$species; length(species1) <- n; species2 <- species$genus; length(species2) <- n; species <- cbind(species1,species2)
strain <- strain$species; length(strain) <- n

taxalevels <- as.data.frame(cbind(strain,species,genus,family,order,class,phylum))
colnames(taxalevels) <- c("strain$species","species$species","species$genus","genus$genus","genus$family","family$family","family$order","order$order","order$class","class$class","class$phylum","phylum$phylum")
for (i in c(1,3,5,7,9,11)){
  ref <- gsub("\\$.*","",colnames(taxalevels)[i+1])
  ref_len <- taxalevels[,i+1]; ref_len <- length(ref_len[!is.na(ref_len)])
  subject <- gsub("\\$.*","",colnames(taxalevels)[i])
  subject_len <- taxalevels[,i]; subject_len <- length(subject_len[!is.na(subject_len)])
  overlap <- intersect(taxalevels[,i+1], taxalevels[,i]); overlap <- overlap[!is.na(overlap)]
  ref_uniq <- taxalevels[,i+1][!c(taxalevels[,i+1]) %in% c(overlap)]; ref_uniq <- length(ref_uniq[!is.na(ref_uniq)])
  subject_uniq <- taxalevels[,i][!c(taxalevels[,i]) %in% c(overlap)]; subject_uniq <- length(subject_uniq[!is.na(subject_uniq)])
  overlap <- length(overlap)
  false_positive <- subject_uniq / subject_len
  false_negative <- ref_uniq / ref_len
  sensitivity <- overlap / (overlap + ref_uniq)
  ctv_hold <- c(ref,ref_len,subject,subject_len,overlap,ref_uniq,subject_uniq,false_positive,false_negative,sensitivity)
  ctv <- rbind(ctv,ctv_hold)
}
colnames(ctv) <- c("ref","ref_len","subject","subject_len","overlap","ref_uniq","subjec_uniq","false_positive","false_negative","sensitivity")


ctv$Cross_taxon_comparison <- paste(ctv$subject,"_vs._",ctv$ref,"_(ref)",sep="") 
ctv$Cross_taxon_comparison <- factor(as.character(ctv$Cross_taxon_comparison), levels=unique(ctv$Cross_taxon_comparison))
write.table(ctv,"cross_taxon_comparison_validated.txt", sep="\t",row.names=FALSE, col.names=T, quote = F)
ctv <- subset(ctv, select=c(Cross_taxon_comparison,false_positive,false_negative,sensitivity))
ctv <- melt(ctv, id="Cross_taxon_comparison")
ctv$value <- as.numeric(ctv$value)

plot <- ggplot(ctv,aes(x=Cross_taxon_comparison,y=value,group=variable,color=variable))+
  geom_line(size=1.5)+ylim(0,1)+
  geom_point()+xlab("Cross Taxon Comparison (validated)")+
  ylab("Value")+ theme_bw() + theme(legend.position="bottom") + 
  theme(text = element_text(size = 20),axis.text.x = element_text(angle=30,size=18.5,hjust = 1)) + labs(color='Metric') +
  theme(plot.margin = margin(2,2,2,2, "cm"))
ggsave(file="cross_taxon_comparison_validated.tiff", plot=plot, width=12, height=6, units=("in"), dpi=300, compression = "lzw")

