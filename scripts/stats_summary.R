#!/usr/bin/env Rscript


args <- commandArgs(TRUE)
sighits <- args[1]
libdir <- args[3]

.libPaths( c( .libPaths(), libdir) )
library(plyr, quietly=T)
library(reshape2, quietly=T)
suppressMessages(library(data.table, quietly=T))




if (args[2] == "strain"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,2,8,9))
    stats1 <- subset(stats1, grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", stats1$abundance))
    stats1$abundance <- as.numeric(as.character(stats1$abundance))
    stats1 <- stats1[!is.na(stats1$abundance),]
  }
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, abundance > "0")
    stats1 <- stats1[!(stats1$sseqid == ""), ]
    stats1 <- stats1[!(stats1$staxids == ""), ]
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$sseqid))
    stats1 <- subset(stats1, grepl('^\\d+$', stats1$staxids))
  }
  if (nrow(stats1) > 0) {
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1$abundance2 <- stats1$abundance
    stats1 <- subset(stats1, select=c(3,1,5))
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats1 <- ddply(stats1, c("staxids"), summarise,
                    N    = length(staxids),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[2] == "species"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,2,13,9))
    stats1 <- subset(stats1, grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", stats1$abundance))
    stats1$abundance <- as.numeric(as.character(stats1$abundance))
    stats1 <- stats1[!is.na(stats1$abundance),]
  }
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, abundance > "0")
    stats1 <- stats1[!(stats1$sseqid == ""), ]
    stats1 <- stats1[!(stats1$species == ""), ]
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$sseqid))
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$species))
  }
  if (nrow(stats1) > 0) {
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1$abundance2 <- stats1$abundance
    stats1 <- subset(stats1, select=c(3,1,5))
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats1 <- ddply(stats1, c("species"), summarise,
                    N    = length(species),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[2] == "genus"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,2,9))
    stats1 <- subset(stats1, grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", stats1$abundance))
    stats1$abundance <- as.numeric(as.character(stats1$abundance))
    stats1 <- stats1[!is.na(stats1$abundance),]
  }
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, abundance > "0")
    stats1 <- stats1[!(stats1$sseqid == ""), ]
    stats1 <- stats1[!(stats1$genus == ""), ]
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$sseqid))
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$genus))
  }
  if (nrow(stats1) > 0) {
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1$abundance2 <- stats1$abundance
    stats1 <- subset(stats1, select=c(3,1,5))
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats1 <- ddply(stats1, c("genus"), summarise,
                    N    = length(genus),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[2] == "family"){
fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,2,13,9))
    stats1 <- subset(stats1, grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", stats1$abundance))
    stats1$abundance <- as.numeric(as.character(stats1$abundance))
    stats1 <- stats1[!is.na(stats1$abundance),]
  }
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, abundance > "0")
    stats1 <- stats1[!(stats1$sseqid == ""), ]
    stats1 <- stats1[!(stats1$family == ""), ]
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$sseqid))
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$family))
  }
  if (nrow(stats1) > 0) {
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1$abundance2 <- stats1$abundance
    stats1 <- subset(stats1, select=c(3,1,5))
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats1 <- ddply(stats1, c("family"), summarise,
                    N    = length(family),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[2] == "order"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,2,13,9))
    stats1 <- subset(stats1, grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", stats1$abundance))
    stats1$abundance <- as.numeric(as.character(stats1$abundance))
    stats1 <- stats1[!is.na(stats1$abundance),]
  }
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, abundance > "0")
    stats1 <- stats1[!(stats1$sseqid == ""), ]
    stats1 <- stats1[!(stats1$order == ""), ]
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$sseqid))
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$order))
  }
  if (nrow(stats1) > 0) {
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1$abundance2 <- stats1$abundance
    stats1 <- subset(stats1, select=c(3,1,5))
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats1 <- ddply(stats1, c("order"), summarise,
                    N    = length(order),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[2] == "class"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,2,13,9))
    stats1 <- subset(stats1, grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", stats1$abundance))
    stats1$abundance <- as.numeric(as.character(stats1$abundance))
    stats1 <- stats1[!is.na(stats1$abundance),]
  }
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, abundance > "0")
    stats1 <- stats1[!(stats1$sseqid == ""), ]
    stats1 <- stats1[!(stats1$class == ""), ]
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$sseqid))
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$class))
  }
  if (nrow(stats1) > 0) {
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1$abundance2 <- stats1$abundance
    stats1 <- subset(stats1, select=c(3,1,5))
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats1 <- ddply(stats1, c("class"), summarise,
                    N    = length(class),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[2] == "phylum"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,2,13,9))
    stats1 <- subset(stats1, grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", stats1$abundance))
    stats1$abundance <- as.numeric(as.character(stats1$abundance))
    stats1 <- stats1[!is.na(stats1$abundance),]
  }
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, abundance > "0")
    stats1 <- stats1[!(stats1$sseqid == ""), ]
    stats1 <- stats1[!(stats1$phylum == ""), ]
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$sseqid))
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$phylum))
  }
  if (nrow(stats1) > 0) {
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1$abundance2 <- stats1$abundance
    stats1 <- subset(stats1, select=c(3,1,5))
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats1 <- ddply(stats1, c("phylum"), summarise,
                    N    = length(phylum),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}
