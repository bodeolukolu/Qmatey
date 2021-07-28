#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

order_mean  <- read.delim(args[1], header = T, sep="\t", check.names=FALSE, fill=TRUE)
mean_agg =aggregate(x=order_mean, by=list(order_mean$order), FUN = "mean")
metag <- subset(mean_agg, select=-c(tax_id,species,genus,family,class,phylum,kingdom,superkingdom))
metag <- subset(metag, select=c(ncol(metag),1:(ncol(metag)-1)))
metag <- subset(metag, select=-c(taxname))
metag <- subset(metag, select=-c(order))
colnames(metag)[colnames(metag)=="Group.1"] <- "order"
write.table(metag, file="order_mean.txt", sep="\t", row.names = FALSE)

order_uniq <- read.delim(args[2], header = T, sep="\t", check.names=FALSE, fill=TRUE)
uniq_agg=aggregate(x=order_uniq, by=list(order_uniq$order), FUN = "mean")
metag_uniq <- subset(uniq_agg, select=-c(tax_id,species,genus,family,class,phylum,kingdom,superkingdom))
metag_uniq <- subset(metag_uniq, select=c(ncol(metag_uniq),1:(ncol(metag_uniq)-1)))
metag_uniq <- subset(metag_uniq, select=-c(taxname))
metag_uniq <- subset(metag_uniq, select=-c(order))
colnames(metag_uniq)[colnames(metag_uniq)=="Group.1"] <- "order"
write.table(metag_uniq, file="order_unique_sequences.txt", sep="\t", row.names = FALSE)

order_quant <- read.delim(args[3], header = T, sep="\t", check.names=FALSE, fill=TRUE)
quant_agg=aggregate(x=order_quant, by=list(order_quant$order), FUN = "mean")
metag_quant <- subset(uniq_agg, select=-c(tax_id,species,genus,family,class,phylum,kingdom,superkingdom))
metag_quant <- subset(metag_quant, select=c(ncol(metag_quant),1:(ncol(metag_quant)-1)))
metag_quant <- subset(metag_quant, select=-c(taxname))
metag_quant <- subset(metag_quant, select=-c(order))
colnames(metag_quant)[colnames(metag_quant)=="Group.1"] <- "family"
write.table(metag_uniq, file="order_quantification_accuracy.txt", sep="\t", row.names = FALSE)