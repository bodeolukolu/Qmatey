#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
ref_taxlevel <- args[3]
libdir <- args[4]
.libPaths( c( .libPaths(), libdir) )
library(dplyr, quietly = T)


ref <- read.delim(file=paste(args[1],"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
ref <- subset(ref, domain != "Viruses")
uniq_seq <- read.delim(file=paste(args[2],"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
uniq_seq <- subset(uniq_seq, domain != "Viruses")


test_uniq <- setdiff(uniq_seq[,colnames(uniq_seq)==ref_taxlevel], ref[,1])
test_uniq <- test_uniq[!is.na(test_uniq)]


uniq_seq <- read.delim(file=paste(args[2],"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
mean <- read.delim(file=paste(args[2],"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
norm_mean <- read.delim(file=paste(args[2],"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
stderr <- read.delim(file=paste(args[2],"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
rel_stderr <- read.delim(file=paste(args[2],"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)

uniq_seq <- uniq_seq[!grepl(paste(test_uniq, collapse="|"), uniq_seq[,colnames(uniq_seq)==ref_taxlevel]),]
mean <- mean[!grepl(paste(test_uniq, collapse="|"), mean[,colnames(mean)==ref_taxlevel]),]
norm_mean <- norm_mean[!grepl(paste(test_uniq, collapse="|"), norm_mean[,colnames(norm_mean)==ref_taxlevel]),]
stderr <- stderr[!grepl(paste(test_uniq, collapse="|"), stderr[,colnames(stderr)==ref_taxlevel]),]
rel_stderr <- rel_stderr[!grepl(paste(test_uniq, collapse="|"), rel_stderr[,colnames(rel_stderr)==ref_taxlevel]),]

write.table(uniq_seq,paste(args[5],"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
write.table(mean,paste(args[5],"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
write.table(norm_mean,paste(args[5],"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
write.table(stderr,paste(args[5],"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
write.table(rel_stderr,paste(args[5],"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)

