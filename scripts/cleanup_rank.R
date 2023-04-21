#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
libdir <- args[2]
.libPaths( c( .libPaths(), libdir) )


uniq_seq <- read.delim(file=paste(args[1],"_taxainfo_unique_sequences.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
mean <- read.delim(file=paste(args[1],"_taxainfo_mean.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
norm_mean <- read.delim(file=paste(args[1],"_taxainfo_mean_normalized.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
stderr <- read.delim(file=paste(args[1],"_taxainfo_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)
rel_stderr <- read.delim(file=paste(args[1],"_taxainfo_rel_quantification_accuracy.txt",sep=""), header=T, sep="\t", fill= T, quote="", check.names = T)


uniq_seq[is.na(uniq_seq)] <- "NA"
mean[is.na(mean)] <- "NA"
norm_mean[is.na(norm_mean)] <- "NA" 
stderr[is.na(stderr)] <- "NA"
rel_stderr[is.na(rel_stderr)] <- "NA"


uniq_seq <- subset(uniq_seq, genus != family)
uniq_seq <- subset(uniq_seq, genus != order)
uniq_seq <- subset(uniq_seq, genus != class)
uniq_seq <- subset(uniq_seq, genus != phylum)


mean <- subset(mean, genus != family)
mean <- subset(mean, genus != order)
mean <- subset(mean, genus != class)
mean <- subset(mean, genus != phylum)
norm_mean <- subset(norm_mean, genus != family)
norm_mean <- subset(norm_mean, genus != order)
norm_mean <- subset(norm_mean, genus != class)
norm_mean <- subset(norm_mean, genus != phylum)
stderr <- subset(stderr, genus != family)
stderr <- subset(stderr, genus != order)
stderr <- subset(stderr, genus != class)
stderr <- subset(stderr, genus != phylum)
rel_stderr <- subset(rel_stderr, genus != family)
rel_stderr <- subset(rel_stderr, genus != order)
rel_stderr <- subset(rel_stderr, genus != class)
rel_stderr <- subset(rel_stderr, genus != phylum)

if(args[1] != "genus") {
  mean <- subset(mean, species != family)
  mean <- subset(mean, species != order)
  mean <- subset(mean, species != class)
  mean <- subset(mean, species != phylum)
  norm_mean <- subset(norm_mean, species != family)
  norm_mean <- subset(norm_mean, species != order)
  norm_mean <- subset(norm_mean, species != class)
  norm_mean <- subset(norm_mean, species != phylum)
  stderr <- subset(stderr, species != family)
  stderr <- subset(stderr, species != order)
  stderr <- subset(stderr, species != class)
  stderr <- subset(stderr, species != phylum)
  rel_stderr <- subset(rel_stderr, species != family)
  rel_stderr <- subset(rel_stderr, species != order)
  rel_stderr <- subset(rel_stderr, species != class)
  rel_stderr <- subset(rel_stderr, species != phylum)
}
  

write.table(uniq_seq,paste(args[1],"_taxainfo_unique_sequences.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
write.table(mean,paste(args[1],"_taxainfo_mean.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
write.table(norm_mean,paste(args[1],"_taxainfo_mean_normalized.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
write.table(stderr,paste(args[1],"_taxainfo_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
write.table(rel_stderr,paste(args[1],"_taxainfo_rel_quantification_accuracy.txt",sep=""), sep="\t",row.names=FALSE, col.names=T, quote = F)
