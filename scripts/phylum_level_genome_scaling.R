#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
taxalevel <- (args[1])
libdir <- args[2]

.libPaths( c( .libPaths(), libdir) )
library(dplyr, quietly = T)

if (taxalevel == "strain"){
  final <- read.delim(file=paste(taxalevel,"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
  final <- subset(final, select=-c(taxname,species,genus,family,order,class,kingdom,domain))
  final$mean_max <- do.call(pmax, final[,2:(ncol(final)-1)])
  final <- subset(final, select=c(tax_id,phylum,mean_max))
  finalq <- final
  finalq[finalq == 0] <- NA
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q99 = quantile(mean_max, probs = 0.98))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q99 * 0.05, digit=0)
  final$threshold <- finalq$threshold[match(final$phylum, finalq$phylum)]
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  keep_tax_id <- final$tax_id
  
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
  final$mean_max <- do.call(pmax, final[,2:(ncol(final)-1)])
  final <- subset(final, select=c(species,phylum,mean_max))
  finalq <- final
  finalq[finalq == 0] <- NA
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q99 = quantile(mean_max, probs = 0.98))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q99 * 0.05, digit=0)
  final$threshold <- finalq$threshold[match(final$phylum, finalq$phylum)]
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  keep_species <- final$species
  
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
  final$mean_max <- do.call(pmax, final[,2:(ncol(final)-1)])
  final <- subset(final, select=c(genus,phylum,mean_max))
  finalq <- final
  finalq[finalq == 0] <- NA
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q99 = quantile(mean_max, probs = 0.98))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q99 * 0.05, digit=0)
  final$threshold <- finalq$threshold[match(final$phylum, finalq$phylum)]
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  keep_genus <- final$genus
  
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
  final$mean_max <- do.call(pmax, final[,2:(ncol(final)-1)])
  final <- subset(final, select=c(family,phylum,mean_max))
  finalq <- final
  finalq[finalq == 0] <- NA
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q99 = quantile(mean_max, probs = 0.98))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q99 * 0.05, digit=0)
  final$threshold <- finalq$threshold[match(final$phylum, finalq$phylum)]
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  keep_family <- final$family
  
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
  final$mean_max <- do.call(pmax, final[,2:(ncol(final)-1)])
  final <- subset(final, select=c(order,phylum,mean_max))
  finalq <- final
  finalq[finalq == 0] <- NA
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q99 = quantile(mean_max, probs = 0.98))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q99 * 0.05, digit=0)
  final$threshold <- finalq$threshold[match(final$phylum, finalq$phylum)]
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  keep_order <- final$order
  
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
  final$mean_max <- do.call(pmax, final[,2:(ncol(final)-1)])
  final <- subset(final, select=c(class,phylum,mean_max))
  finalq <- final
  finalq[finalq == 0] <- NA
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q99 = quantile(mean_max, probs = 0.98))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q99 * 0.05, digit=0)
  final$threshold <- finalq$threshold[match(final$phylum, finalq$phylum)]
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  keep_class <- final$class
  
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
  final$mean_max <- do.call(pmax, final[,2:ncol(final)])
  final <- subset(final, select=c(phylum,mean_max))
  finalq <- final
  finalq[finalq == 0] <- NA
  finalq <- finalq %>%
    group_by(phylum) %>%
    summarize(q99 = quantile(mean_max, probs = 0.98))
  finalq <- as.data.frame(finalq)
  finalq$threshold <- round(finalq$q99 * 0.05, digit=0)
  final$threshold <- finalq$threshold[match(final$phylum, finalq$phylum)]
  final$keep <- final$mean_max - final$threshold
  final <- subset(final, final$keep > 0)
  keep_phylum <- final$phylum
  
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


