#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

class_mean  <- read.delim(args[1], header = T, sep="\t", check.names=FALSE, fill=TRUE)
mean_agg =aggregate(x=class_mean, by=list(class_mean$class), FUN = "mean")
metag <- subset(mean_agg, select=-c(tax_id,species,genus,family,order,phylum,kingdom,superkingdom))
metag <- subset(metag, select=c(ncol(metag),1:(ncol(metag)-1)))
metag <- subset(metag, select=-c(taxname))
metag <- subset(metag, select=-c(class))
colnames(metag)[colnames(metag)=="Group.1"] <- "class"
write.table(metag, file="class_mean.txt", sep="\t", row.names = FALSE)

class_uniq <- read.delim(args[2], header = T, sep="\t", check.names=FALSE, fill=TRUE)
uniq_agg=aggregate(x=class_uniq, by=list(class_uniq$class), FUN = "mean")
metag_uniq <- subset(uniq_agg, select=-c(tax_id,species,genus,family,order,phylum,kingdom,superkingdom))
metag_uniq <- subset(metag_uniq, select=c(ncol(metag_uniq),1:(ncol(metag_uniq)-1)))
metag_uniq <- subset(metag_uniq, select=-c(taxname))
metag_uniq <- subset(metag_uniq, select=-c(class))
colnames(metag_uniq)[colnames(metag_uniq)=="Group.1"] <- "class"
write.table(metag_uniq, file="class_unique_sequences.txt", sep="\t", row.names = FALSE)

class_quant <- read.delim(args[3], header = T, sep="\t", check.names=FALSE, fill=TRUE)
quant_agg=aggregate(x=class_quant, by=list(class_quant$class), FUN = "mean")
metag_quant <- subset(uniq_agg, select=-c(tax_id,species,genus,family,order,phylum,kingdom,superkingdom))
metag_quant <- subset(metag_quant, select=c(ncol(metag_quant),1:(ncol(metag_quant)-1)))
metag_quant <- subset(metag_quant, select=-c(taxname))
metag_quant <- subset(metag_quant, select=-c(class))
colnames(metag_quant)[colnames(metag_quant)=="Group.1"] <- "class"
write.table(metag_uniq, file="class_quantification_accuracy.txt", sep="\t", row.names = FALSE)