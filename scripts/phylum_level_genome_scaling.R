#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
taxalevel <- (args[1])
libdir <- args[2]
strain_min_uniq_thresh <- as.numeric(args[3])
zero_inflated <- as.numeric(args[4])
wgs <- args[5]
.libPaths( c( .libPaths(), libdir) )
library(dplyr, quietly = T)

if (taxalevel == "strain"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(taxname,species,genus,family,order,class,kingdom,domain))
  final$mean_mean <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, mean)
  final <- subset(final, select=c(tax_id,phylum,mean_mean))
  finalq <- final
  finalq$mean_mean[finalq$mean_mean < strain_min_uniq_thresh] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_mean, probs = 0.95), q5 = quantile(mean_mean, probs = 0.05))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_mean - final$threshold
  # final <- subset(final, final$keep >= 0)
  if (wgs == "true") { final <- subset(final, final$mean_mean > final$q5) }
  keep_tax_id1 <- final$tax_id
  
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- t(subset(final, select=-c(tax_id,taxname,species,genus,family,order,class,phylum,kingdom,domain)))
  final[final <= strain_min_uniq_thresh] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final))
  final <-  subset(final, final$zinflate <= zero_inflated)
  rm_samples <- row.names(final)
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(taxname,species,genus,family,order,class,kingdom,domain))
  final[,1:(ncol(final)-1)] [final[,1:(ncol(final)-1)] <= 1] <- NA
  final$miss <- rowSums(is.na(final))
  sampleN <- ncol(final)-3
  final <- subset(final, select=c(tax_id,phylum,miss))
  final <- subset(final, final$miss < sampleN)
  keep_tax_id2 <- final$tax_id
  keep_tax_id <- intersect(keep_tax_id1, keep_tax_id2)
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$tax_id %in% keep_tax_id,]
  uniq_seq <- uniq_seq[ !rownames(uniq_seq) %in% c(rm_samples),]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$tax_id %in% keep_tax_id,]
  mean <- mean[ !rownames(mean) %in% c(rm_samples),]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$tax_id %in% keep_tax_id,]
  norm_mean <- norm_mean[ !rownames(norm_mean) %in% c(rm_samples),]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$tax_id %in% keep_tax_id,]
  stderr <- stderr[ rownames(stderr) %in% c(rm_samples),]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$tax_id %in% keep_tax_id,]
  rel_stderr <- rel_stderr[ !rownames(rel_stderr) %in% c(rm_samples),]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "species"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(genus,family,order,class,kingdom,domain))
  final$mean_mean <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, mean)
  final <- subset(final, select=c(species,phylum,mean_mean))
  finalq <- final
  finalq$mean_mean[finalq$mean_mean < 2] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_mean, probs = 0.95), q5 = quantile(mean_mean, probs = 0.05))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_mean - final$threshold
  # final <- subset(final, final$keep > 0)
  if (wgs == "true") { final <- subset(final, final$mean_mean > final$q5) }
  keep_species1 <- final$species
  
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- t(subset(final, select=-c(species,genus,family,order,class,phylum,kingdom,domain)))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final))
  final <-  subset(final, final$zinflate <= zero_inflated)
  rm_samples <- row.names(final)
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(genus,family,order,class,kingdom,domain))
  final[,1:(ncol(final)-1)] [final[,1:(ncol(final)-1)] <= 1] <- NA
  final$miss <- rowSums(is.na(final))
  sampleN <- ncol(final)-3
  final <- subset(final, select=c(species,phylum,miss))
  final <- subset(final, final$miss < sampleN)
  keep_species2 <- final$species
  keep_species <- intersect(keep_species1, keep_species2)
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$species %in% keep_species,]
  uniq_seq <- uniq_seq[ !rownames(uniq_seq) %in% c(rm_samples),]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$species %in% keep_species,]
  mean <- mean[ !rownames(mean) %in% c(rm_samples),]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$species %in% keep_species,]
  norm_mean <- norm_mean[ !rownames(norm_mean) %in% c(rm_samples),]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$species %in% keep_species,]
  stderr <- stderr[ rownames(stderr) %in% c(rm_samples),]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$species %in% keep_species,]
  rel_stderr <- rel_stderr[ !rownames(rel_stderr) %in% c(rm_samples),]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "genus"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(family,order,class,kingdom,domain))
  final$mean_mean <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, mean)
  final <- subset(final, select=c(genus,phylum,mean_mean))
  finalq <- final
  finalq$mean_mean[finalq$mean_mean < 2] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_mean, probs = 0.95), q5 = quantile(mean_mean, probs = 0.05))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_mean - final$threshold
  # final <- subset(final, final$keep >= 0)
  if (wgs == "true") { final <- subset(final, final$mean_mean > final$q5) }
  keep_genus1 <- final$genus
  
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- t(subset(final, select=-c(genus,family,order,class,phylum,kingdom,domain)))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final))
  final <-  subset(final, final$zinflate <= zero_inflated)
  rm_samples <- row.names(final)
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(family,order,class,kingdom,domain))
  final[,1:(ncol(final)-1)] [final[,1:(ncol(final)-1)] <= 1] <- NA
  final$miss <- rowSums(is.na(final))
  sampleN <- ncol(final)-3
  final <- subset(final, select=c(genus,phylum,miss))
  final <- subset(final, final$miss < sampleN)
  keep_genus2 <- final$genus
  keep_genus <- intersect(keep_genus1, keep_genus2)
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$genus %in% keep_genus,]
  uniq_seq <- uniq_seq[ !rownames(uniq_seq) %in% c(rm_samples),]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$genus %in% keep_genus,]
  mean <- mean[ !rownames(mean) %in% c(rm_samples),]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$genus %in% keep_genus,]
  norm_mean <- norm_mean[ !rownames(norm_mean) %in% c(rm_samples),]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$genus %in% keep_genus,]
  stderr <- stderr[ rownames(stderr) %in% c(rm_samples),]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$genus %in% keep_genus,]
  rel_stderr <- rel_stderr[ !rownames(rel_stderr) %in% c(rm_samples),]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "family"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(order,class,kingdom,domain))
  final$mean_mean <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, mean)
  final <- subset(final, select=c(family,phylum,mean_mean))
  finalq <- final
  finalq$mean_mean[finalq$mean_mean < 2] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_mean, probs = 0.95), q5 = quantile(mean_mean, probs = 0.05))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_mean - final$threshold
  # final <- subset(final, final$keep >= 0)
  if (wgs == "true") { final <- subset(final, final$mean_mean > final$q5) }
  keep_family1 <- final$family
  
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- t(subset(final, select=-c(family,order,class,phylum,kingdom,domain)))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final))
  final <-  subset(final, final$zinflate <= zero_inflated)
  rm_samples <- row.names(final)
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(order,class,kingdom,domain))
  final[,1:(ncol(final)-1)] [final[,1:(ncol(final)-1)] <= 1] <- NA
  final$miss <- rowSums(is.na(final))
  sampleN <- ncol(final)-3
  final <- subset(final, select=c(family,phylum,miss))
  final <- subset(final, final$miss < sampleN)
  keep_family2 <- final$family
  keep_family <- intersect(keep_family1, keep_family2)
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$family %in% keep_family,]
  uniq_seq <- uniq_seq[ !rownames(uniq_seq) %in% c(rm_samples),]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$family %in% keep_family,]
  mean <- mean[ !rownames(mean) %in% c(rm_samples),]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$family %in% keep_family,]
  norm_mean <- norm_mean[ !rownames(norm_mean) %in% c(rm_samples),]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$family %in% keep_family,]
  stderr <- stderr[ rownames(stderr) %in% c(rm_samples),]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$family %in% keep_family,]
  rel_stderr <- rel_stderr[ !rownames(rel_stderr) %in% c(rm_samples),]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "order"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(class,kingdom,domain))
  final$mean_mean <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, mean)
  final <- subset(final, select=c(order,phylum,mean_mean))
  finalq <- final
  finalq$mean_mean[finalq$mean_mean < 2] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_mean, probs = 0.95), q5 = quantile(mean_mean, probs = 0.05))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_mean - final$threshold
  # final <- subset(final, final$keep >= 0)
  if (wgs == "true") { final <- subset(final, final$mean_mean > final$q5) }
  keep_order1 <- final$order
  
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- t(subset(final, select=-c(order,class,phylum,kingdom,domain)))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final))
  final <-  subset(final, final$zinflate <= zero_inflated)
  rm_samples <- row.names(final)
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(class,kingdom,domain))
  final[,1:(ncol(final)-1)] [final[,1:(ncol(final)-1)] <= 1] <- NA
  final$miss <- rowSums(is.na(final))
  sampleN <- ncol(final)-3
  final <- subset(final, select=c(order,phylum,miss))
  final <- subset(final, final$miss < sampleN)
  keep_order2 <- final$order
  keep_order <- intersect(keep_order1, keep_order2)
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$order %in% keep_order,]
  uniq_seq <- uniq_seq[ !rownames(uniq_seq) %in% c(rm_samples),]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$order %in% keep_order,]
  mean <- mean[ !rownames(mean) %in% c(rm_samples),]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$order %in% keep_order,]
  norm_mean <- norm_mean[ !rownames(norm_mean) %in% c(rm_samples),]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$order %in% keep_order,]
  stderr <- stderr[ rownames(stderr) %in% c(rm_samples),]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$order %in% keep_order,]
  rel_stderr <- rel_stderr[ !rownames(rel_stderr) %in% c(rm_samples),]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "class"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(kingdom,domain))
  final$mean_mean <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, mean)
  final <- subset(final, select=c(class,phylum,mean_mean))
  finalq <- final
  finalq$mean_mean[finalq$mean_mean < 2] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_mean, probs = 0.95), q5 = quantile(mean_mean, probs = 0.05))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_mean - final$threshold
  # final <- subset(final, final$keep >= 0)
  if (wgs == "true") { final <- subset(final, final$mean_mean > final$q5) }
  keep_class1 <- final$class
  
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- (subset(final, select=-c(class,phylum,kingdom,domain)))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final))
  final <- subset(final, final$zinflate <= zero_inflated)
  rm_samples <- row.names(final)
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(kingdom,domain))
  final[,1:(ncol(final)-1)] [final[,1:(ncol(final)-1)] <= 1] <- NA
  final$miss <- rowSums(is.na(final))
  sampleN <- ncol(final)-3
  final <- subset(final, select=c(class,phylum,miss))
  final <- subset(final, final$miss < sampleN)
  keep_class2 <- final$class
  keep_class <- intersect(keep_class1, keep_class2)
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$class %in% keep_class,]
  uniq_seq <- uniq_seq[!rownames(uniq_seq) %in% c(rm_samples),]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$class %in% keep_class,]
  mean <- mean[ !rownames(mean) %in% c(rm_samples),]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$class %in% keep_class,]
  norm_mean <- norm_mean[ !rownames(norm_mean) %in% c(rm_samples),]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$class %in% keep_class,]
  stderr <- stderr[ rownames(stderr) %in% c(rm_samples),]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$class %in% keep_class,]
  rel_stderr <- rel_stderr[ !rownames(rel_stderr) %in% c(rm_samples),]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

if (taxalevel == "phylum"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(phylum,kingdom,domain))
  final[final <= 2] <- NA; final <- as.data.frame(final)
  final$zinflate <- rowSums(is.na(final))/(ncol(final))
  final <- subset(final, final$zinflate <= zero_inflated)
  rm_samples <- row.names(final)
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(kingdom,domain))
  final[,1:ncol(final)] [final[,1:ncol(final)] <= 1] <- NA
  final$miss <- rowSums(is.na(final))
  sampleN <- ncol(final)-2
  final <- subset(final, select=c(phylum,miss))
  final <- subset(final, final$miss < sampleN)
  keep_phylum <- final$phylum
  
  
  uniq_seq <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  norm_mean <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  stderr <- read.delim(file=paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  rel_stderr <- read.delim(file=paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  
  uniq_seq <- uniq_seq[uniq_seq$phylum %in% keep_phylum,]
  uniq_seq <- uniq_seq[!rownames(uniq_seq) %in% c(rm_samples),]
  write.table(uniq_seq,paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  mean <- mean[mean$phylum %in% keep_phylum,]
  mean <- mean[!rownames(mean) %in% c(rm_samples),]
  write.table(mean,paste(taxalevel,"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  norm_mean <- norm_mean[norm_mean$phylum %in% keep_phylum,]
  norm_mean <- norm_mean[!rownames(norm_mean) %in% c(rm_samples),]
  write.table(norm_mean,paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  stderr <- stderr[stderr$phylum %in% keep_phylum,]
  stderr <- stderr[!rownames(stderr) %in% c(rm_samples),]
  write.table(stderr,paste(taxalevel,"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
  rel_stderr <- rel_stderr[rel_stderr$phylum %in% keep_phylum,]
  rel_stderr <- rel_stderr[!rownames(rel_stderr) %in% c(rm_samples),]
  write.table(rel_stderr,paste(taxalevel,"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
}

