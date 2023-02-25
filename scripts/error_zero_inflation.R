#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
taxalevel <- (args[1])
libdir <- args[2]
strain_min_uniq_thresh <- as.numeric(args[3])
zero_inflated <- as.numeric(args[4])
.libPaths( c( .libPaths(), libdir) )
library(dplyr, quietly = T)

if (taxalevel == "strain"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(taxname,species,genus,family,order,class,phylum,kingdom,domain))
  final[final <= strain_min_uniq_thresh] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final)-1)
  final <-  subset(final, final$zinflate > zero_inflated)
  keep_tax_id <- final$zinflate
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$tax_id %in% keep_tax_id,]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$tax_id %in% keep_tax_id,]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$tax_id %in% keep_tax_id,]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$tax_id %in% keep_tax_id,]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$tax_id %in% keep_tax_id,]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "species"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- t(subset(final, select=-c(genus,family,order,class,phylum,kingdom,domain)))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final)-1)
  final <-  subset(final, final$zinflate > zero_inflated)
  keep_species <- final$zinflate
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$species %in% keep_species,]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$species %in% keep_species,]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$species %in% keep_species,]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$species %in% keep_species,]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$species %in% keep_species,]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "genus"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(family,order,class,phylum,kingdom,domain))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final)-1)
  final <-  subset(final, final$zinflate > zero_inflated)
  keep_genus <- final$zinflate
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$genus %in% keep_genus,]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$genus %in% keep_genus,]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$genus %in% keep_genus,]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$genus %in% keep_genus,]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$genus %in% keep_genus,]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "family"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(order,class,phylum,kingdom,domain))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final)-1)
  final <-  subset(final, final$zinflate > zero_inflated)
  keep_family <- final$zinflate
  
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$family %in% keep_family,]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$family %in% keep_family,]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$family %in% keep_family,]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$family %in% keep_family,]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$family %in% keep_family,]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "order"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(class,phylum,kingdom,domain))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final)-1)
  final <-  subset(final, final$zinflate > zero_inflated)
  keep_order <- final$zinflate  

  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$order %in% keep_order,]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$order %in% keep_order,]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$order %in% keep_order,]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$order %in% keep_order,]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$order %in% keep_order,]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "class"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(phylum,kingdom,domain))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final)-1)
  final <-  subset(final, final$zinflate > zero_inflated)
  keep_class <- final$zinflate  

  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$class %in% keep_class,]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$class %in% keep_class,]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$class %in% keep_class,]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$class %in% keep_class,]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$class %in% keep_class,]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "phylum"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(kingdom,domain))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final[,1:ncol(final)] [final[,1:ncol(final)] <= 1] <- NA
  final$zinflate <- rowSums(is.na(final))/(ncol(final)-1)
  final <-  subset(final, final$zinflate > zero_inflated)
  keep_phylum <- final$zinflate  
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$phylum %in% keep_phylum,]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$phylum %in% keep_phylum,]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$phylum %in% keep_phylum,]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$phylum %in% keep_phylum,]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$phylum %in% keep_phylum,]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

