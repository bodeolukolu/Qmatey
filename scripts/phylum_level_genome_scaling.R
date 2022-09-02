#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
taxalevel <- (args[1])
libdir <- args[2]
strain_min_uniq_thresh <- as.numeric(args[3])
.libPaths( c( .libPaths(), libdir) )
library(dplyr, quietly = T)

if (taxalevel == "strain"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(taxname,species,genus,family,order,class,kingdom,domain))
  final$mean_max <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, max)
  final <- subset(final, select=c(tax_id,phylum,mean_max))
  finalq <- final
  finalq[finalq < strain_min_uniq_thresh] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_max, probs = 0.95), q5 = quantile(mean_max, probs = 0.05))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  final <- subset(final, final$mean_max >= final$q5)
  keep_tax_id1 <- final$tax_id
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
    final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(taxname,species,genus,family,order,class,kingdom,domain))
  final[,2:(ncol(final)-1)] [final[,2:(ncol(final)-1)] <= 1] <- NA
  final$mean_max <- rowSums(is.na(final))
  sampleN <- ncol(final)-2
  final <- subset(final, select=c(tax_id,mean_max))
  final <- subset(final, final$mean_max > 0)
  final <- subset(final, final$mean_max > (sampleN/100))
  keep_tax_id2 <- final$tax_id
  keep_tax_id <- intersect(keep_tax_id1, keep_tax_id2)
  
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
  final <- subset(final, select=-c(genus,family,order,class,kingdom,domain))
  final$mean_max <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, max)
  final <- subset(final, select=c(species,phylum,mean_max))
  finalq <- final
  finalq[finalq < 2] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_max, probs = 0.95), q5 = quantile(mean_max, probs = 0.1))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  final <- subset(final, final$mean_max >= final$q5)
  keep_species1 <- final$species
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(genus,family,order,class,kingdom,domain))
  final[,2:(ncol(final)-1)] [final[,2:(ncol(final)-1)] <= 1] <- NA
  final$mean_max <- rowSums(is.na(final))
  sampleN <- ncol(final)-2
  final <- subset(final, select=c(tax_id,mean_max))
  final <- subset(final, final$mean_max > 0)
  final <- subset(final, final$mean_max > (sampleN/100))
  keep_species2 <- final$species
  keep_species <- intersect(keep_species1, keep_species2)
  
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
  final <- subset(final, select=-c(family,order,class,kingdom,domain))
  final$mean_max <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, max)
  final <- subset(final, select=c(genus,phylum,mean_max))
  finalq <- final
  finalq[finalq < 3] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_max, probs = 0.95), q5 = quantile(mean_max, probs = 0.15))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  final <- subset(final, final$mean_max >= final$q5)
  keep_genus1 <- final$genus
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(family,order,class,kingdom,domain))
  final[,2:(ncol(final)-1)] [final[,2:(ncol(final)-1)] <= 1] <- NA
  final$mean_max <- rowSums(is.na(final))
  sampleN <- ncol(final)-2
  final <- subset(final, select=c(tax_id,mean_max))
  final <- subset(final, final$mean_max > 0)
  final <- subset(final, final$mean_max > (sampleN/100))
  keep_genus2 <- final$genus
  keep_genus <- intersect(keep_genus1, keep_genus2)
  
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
  final <- subset(final, select=-c(order,class,kingdom,domain))
  final$mean_max <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, max)
  final <- subset(final, select=c(family,phylum,mean_max))
  finalq <- final
  finalq[finalq < 4] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_max, probs = 0.95), q5 = quantile(mean_max, probs = 0.2))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  final <- subset(final, final$mean_max >= final$q5)
  keep_family1 <- final$family
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(order,class,kingdom,domain))
  final[,2:(ncol(final)-1)] [final[,2:(ncol(final)-1)] <= 1] <- NA
  final$mean_max <- rowSums(is.na(final))
  sampleN <- ncol(final)-2
  final <- subset(final, select=c(tax_id,mean_max))
  final <- subset(final, final$mean_max > 0)
  final <- subset(final, final$mean_max > (sampleN/100))
  keep_family2 <- final$family
  keep_family <- intersect(keep_family1, keep_family2)
  
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
  final <- subset(final, select=-c(class,kingdom,domain))
  final$mean_max <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, max)
  final <- subset(final, select=c(order,phylum,mean_max))
  finalq <- final
  finalq[finalq < 5] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_max, probs = 0.95), q5 = quantile(mean_max, probs = 0.25))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  final <- subset(final, final$mean_max >= final$q5)
  keep_order1 <- final$order
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(class,kingdom,domain))
  final[,2:(ncol(final)-1)] [final[,2:(ncol(final)-1)] <= 1] <- NA
  final$mean_max <- rowSums(is.na(final))
  sampleN <- ncol(final)-2
  final <- subset(final, select=c(tax_id,mean_max))
  final <- subset(final, final$mean_max > 0)
  final <- subset(final, final$mean_max > (sampleN/100))
  keep_order2 <- final$order
  keep_order <- intersect(keep_order1, keep_order2)
  
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
  final <- subset(final, select=-c(kingdom,domain))
  final$mean_max <- apply(as.matrix(final[,2:(ncol(final)-1)]), 1, max)
  final <- subset(final, select=c(class,phylum,mean_max))
  finalq <- final
  finalq[finalq < 6] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_max, probs = 0.95), q5 = quantile(mean_max, probs = 0.3))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  final <- subset(final, final$mean_max >= final$q5)
  keep_class1 <- final$class
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(kingdom,domain))
  final[,2:(ncol(final)-1)] [final[,2:(ncol(final)-1)] <= 1] <- NA
  final$mean_max <- rowSums(is.na(final))
  sampleN <- ncol(final)-2
  final <- subset(final, select=c(tax_id,mean_max))
  final <- subset(final, final$mean_max > 0)
  final <- subset(final, final$mean_max > (sampleN/100))
  keep_class2 <- final$class
  keep_class <- intersect(keep_class1, keep_class2)
  
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
  final$mean_max <- apply(as.matrix(final[,2:(ncol(final))]), 1, max)
  final <- subset(final, select=c(phylum,mean_max))
  finalq <- final
  finalq[finalq < 7] <- NA
  finalq <- na.omit(finalq)
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q95 = quantile(mean_max, probs = 0.95))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q95 * 0.05, digit=0)
  final <- merge(final, finalq, by=c("phylum"), all.y = T)
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  keep_phylum1 <- final$phylum
  
  if (file.exists(paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""))) {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  } else {
      final <- read.delim(file=paste(taxalevel,"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  }
  final <- subset(final, select=-c(kingdom,domain))
  final[,2:ncol(final)] [final[,2:ncol(final)] <= 1] <- NA
  final$mean_max <- rowSums(is.na(final))
  sampleN <- ncol(final)-1
  final <- subset(final, select=c(tax_id,mean_max))
  final <- subset(final, final$mean_max > 0)
  final <- subset(final, final$mean_max > (sampleN/100))
  keep_phylum2 <- final$phylum
  keep_phylum <- intersect(keep_phylum1, keep_phylum2)
  
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


