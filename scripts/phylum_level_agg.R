#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

phylum_mean  <- read.delim(args[1], header = T, sep="\t", check.names=FALSE, fill=TRUE)
mean_agg =aggregate(x=phylum_mean, by=list(phylum_mean$phylum), FUN = "mean")
metag <- subset(mean_agg, select=-c(tax_id,species,genus,family,order,class,kingdom,superkingdom))
metag <- subset(metag, select=c(ncol(metag),1:(ncol(metag)-1)))
metag <- subset(metag, select=-c(taxname))
metag <- subset(metag, select=-c(phylum))
colnames(metag)[colnames(metag)=="Group.1"] <- "phylum"
write.table(metag, file="phylum_mean.txt", sep="\t", row.names = FALSE)

phylum_uniq <- read.delim(args[2], header = T, sep="\t", check.names=FALSE, fill=TRUE)
uniq_agg=aggregate(x=phylum_uniq, by=list(phylum_uniq$phylum), FUN = "mean")
metag_uniq <- subset(uniq_agg, select=-c(tax_id,species,genus,family,order,class,kingdom,superkingdom))
metag_uniq <- subset(metag_uniq, select=c(ncol(metag_uniq),1:(ncol(metag_uniq)-1)))
metag_uniq <- subset(metag_uniq, select=-c(taxname))
metag_uniq <- subset(metag_uniq, select=-c(phylum))
colnames(metag_uniq)[colnames(metag_uniq)=="Group.1"] <- "phylum"
write.table(metag_uniq, file="phylum_unique_sequences.txt", sep="\t", row.names = FALSE)

phylum_quant <- read.delim(args[3], header = T, sep="\t", check.names=FALSE, fill=TRUE)
quant_agg=aggregate(x=phylum_quant, by=list(phylum_quant$class), FUN = "mean")
metag_quant <- subset(uniq_agg, select=-c(tax_id,species,genus,family,order,class,kingdom,superkingdom))
metag_quant <- subset(metag_quant, select=c(ncol(metag_quant),1:(ncol(metag_quant)-1)))
metag_quant <- subset(metag_quant, select=-c(taxname))
metag_quant <- subset(metag_quant, select=-c(phylum))
colnames(metag_quant)[colnames(metag_quant)=="Group.1"] <- "phylum"
write.table(metag_uniq, file="phylum_quantification_accuracy.txt", sep="\t", row.names = FALSE)