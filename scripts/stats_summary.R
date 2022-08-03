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
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,stitle=stats1$stitle,staxids=stats1$staxids), FUN=sum)
    colnames(stats1)[4] <- c("abundance")
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance2 <- stats1$abundance
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1 <- subset(stats1, select=c(3,4,5))
    stats1$uniqseq <- 0
    stats1$uniqseq[stats1$abundance >= 0] <- 1
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats2 <- ddply(stats1, c("staxids","uniqseq"), summarise, Nt = length(staxids))
    stats3 <- as.data.frame(table(stats2$staxids)); colnames(stats3)[1] <- c("staxids")
    stats1 <- merge(stats1, stats3, by=c("staxids"), all.x=TRUE)
    stats1$abundance[stats1$uniqseq == "0" & stats1$Freq == "2" ] <- NA
    stats1 <- ddply(stats1, c("staxids"), summarise,
                    N    = length(staxids),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE)
    )
    stats2 <- merge(stats2, stats3, by=c("staxids"), all.x=TRUE)
    stats2$Nt[stats2$uniqseq == "0" & stats2$Freq == "2" ] <- NA; stats2 <- na.omit(stats2)
    stats1 <- merge(stats1, stats2, by=c("staxids"), all.x=TRUE)
    stats1$se <- stats1$sd/stats1$Nt    
    stats1 <- subset(stats1, select=c(1,3,2,8))
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
  }
  if (nrow(stats1) > 0) {
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,stitle=stats1$stitle,species=stats1$species), FUN=sum)
    colnames(stats1)[4] <- c("abundance")
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance2 <- stats1$abundance
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1 <- subset(stats1, select=c(3,4,5))
    stats1$uniqseq <- 0
    stats1$uniqseq[stats1$abundance >= 0] <- 1
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats2 <- ddply(stats1, c("species","uniqseq"), summarise, Nt = length(species))
    stats3 <- as.data.frame(table(stats2$species)); colnames(stats3)[1] <- c("species")
    stats1 <- merge(stats1, stats3, by=c("species"), all.x=TRUE)
    stats1$abundance[stats1$uniqseq == "0" & stats1$Freq == "2" ] <- NA
    stats1 <- ddply(stats1, c("species"), summarise,
                    N    = length(species),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE)
    )
    stats2 <- merge(stats2, stats3, by=c("species"), all.x=TRUE)
    stats2$Nt[stats2$uniqseq == "0" & stats2$Freq == "2" ] <- NA; stats2 <- na.omit(stats2)
    stats1 <- merge(stats1, stats2, by=c("species"), all.x=TRUE)
    stats1$se <- stats1$sd/stats1$Nt    
    stats1 <- subset(stats1, select=c(1,3,2,8))
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
    stats1 <- subset(stats1, select=c(1,2,13,9))
    stats1 <- subset(stats1, grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", stats1$abundance))
    stats1$abundance <- as.numeric(as.character(stats1$abundance))
    stats1 <- stats1[!is.na(stats1$abundance),]
  }
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, abundance > "0")
    stats1 <- stats1[!(stats1$sseqid == ""), ]
    stats1 <- stats1[!(stats1$genus == ""), ]
    stats1 <- subset(stats1, !grepl('^\\d+$', stats1$sseqid))
  }
  if (nrow(stats1) > 0) {
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,stitle=stats1$stitle,genus=stats1$genus), FUN=sum)
    colnames(stats1)[4] <- c("abundance")
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance2 <- stats1$abundance
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1 <- subset(stats1, select=c(3,4,5))
    stats1$uniqseq <- 0
    stats1$uniqseq[stats1$abundance >= 0] <- 1
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats2 <- ddply(stats1, c("genus","uniqseq"), summarise, Nt = length(genus))
    stats3 <- as.data.frame(table(stats2$genus)); colnames(stats3)[1] <- c("genus")
    stats1 <- merge(stats1, stats3, by=c("genus"), all.x=TRUE)
    stats1$abundance[stats1$uniqseq == "0" & stats1$Freq == "2" ] <- NA
    stats1 <- ddply(stats1, c("genus"), summarise,
                    N    = length(genus),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE)
    )
    stats2 <- merge(stats2, stats3, by=c("genus"), all.x=TRUE)
    stats2$Nt[stats2$uniqseq == "0" & stats2$Freq == "2" ] <- NA; stats2 <- na.omit(stats2)
    stats1 <- merge(stats1, stats2, by=c("genus"), all.x=TRUE)
    stats1$se <- stats1$sd/stats1$Nt    
    stats1 <- subset(stats1, select=c(1,3,2,8))
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
  }
  if (nrow(stats1) > 0) {
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,stitle=stats1$stitle,family=stats1$family), FUN=sum)
    colnames(stats1)[4] <- c("abundance")
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance2 <- stats1$abundance
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1 <- subset(stats1, select=c(3,4,5))
    stats1$uniqseq <- 0
    stats1$uniqseq[stats1$abundance >= 0] <- 1
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats2 <- ddply(stats1, c("family","uniqseq"), summarise, Nt = length(family))
    stats3 <- as.data.frame(table(stats2$family)); colnames(stats3)[1] <- c("family")
    stats1 <- merge(stats1, stats3, by=c("family"), all.x=TRUE)
    stats1$abundance[stats1$uniqseq == "0" & stats1$Freq == "2" ] <- NA
    stats1 <- ddply(stats1, c("family"), summarise,
                    N    = length(family),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE)
    )
    stats2 <- merge(stats2, stats3, by=c("family"), all.x=TRUE)
    stats2$Nt[stats2$uniqseq == "0" & stats2$Freq == "2" ] <- NA; stats2 <- na.omit(stats2)
    stats1 <- merge(stats1, stats2, by=c("family"), all.x=TRUE)
    stats1$se <- stats1$sd/stats1$Nt    
    stats1 <- subset(stats1, select=c(1,3,2,8))
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
  }
  if (nrow(stats1) > 0) {
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,stitle=stats1$stitle,order=stats1$order), FUN=sum)
    colnames(stats1)[4] <- c("abundance")
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance2 <- stats1$abundance
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1 <- subset(stats1, select=c(3,4,5))
    stats1$uniqseq <- 0
    stats1$uniqseq[stats1$abundance >= 0] <- 1
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats2 <- ddply(stats1, c("order","uniqseq"), summarise, Nt = length(order))
    stats3 <- as.data.frame(table(stats2$order)); colnames(stats3)[1] <- c("order")
    stats1 <- merge(stats1, stats3, by=c("order"), all.x=TRUE)
    stats1$abundance[stats1$uniqseq == "0" & stats1$Freq == "2" ] <- NA
    stats1 <- ddply(stats1, c("order"), summarise,
                    N    = length(order),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE)
    )
    stats2 <- merge(stats2, stats3, by=c("order"), all.x=TRUE)
    stats2$Nt[stats2$uniqseq == "0" & stats2$Freq == "2" ] <- NA; stats2 <- na.omit(stats2)
    stats1 <- merge(stats1, stats2, by=c("order"), all.x=TRUE)
    stats1$se <- stats1$sd/stats1$Nt    
    stats1 <- subset(stats1, select=c(1,3,2,8))
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
  }
  if (nrow(stats1) > 0) {
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,stitle=stats1$stitle,class=stats1$class), FUN=sum)
    colnames(stats1)[4] <- c("abundance")
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance2 <- stats1$abundance
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1 <- subset(stats1, select=c(3,4,5))
    stats1$uniqseq <- 0
    stats1$uniqseq[stats1$abundance >= 0] <- 1
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats2 <- ddply(stats1, c("class","uniqseq"), summarise, Nt = length(class))
    stats3 <- as.data.frame(table(stats2$class)); colnames(stats3)[1] <- c("class")
    stats1 <- merge(stats1, stats3, by=c("class"), all.x=TRUE)
    stats1$abundance[stats1$uniqseq == "0" & stats1$Freq == "2" ] <- NA
    stats1 <- ddply(stats1, c("class"), summarise,
                    N    = length(class),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE)
    )
    stats2 <- merge(stats2, stats3, by=c("class"), all.x=TRUE)
    stats2$Nt[stats2$uniqseq == "0" & stats2$Freq == "2" ] <- NA; stats2 <- na.omit(stats2)
    stats1 <- merge(stats1, stats2, by=c("class"), all.x=TRUE)
    stats1$se <- stats1$sd/stats1$Nt    
    stats1 <- subset(stats1, select=c(1,3,2,8))
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
  }
  if (nrow(stats1) > 0) {
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,stitle=stats1$stitle,phylum=stats1$phylum), FUN=sum)
    colnames(stats1)[4] <- c("abundance")
    stats1$stitle <- gsub(".*(ribosomal_RNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(ribosomal_RNA).*", "\\1", stats1$stitle)
    stats1$stitle <- gsub(".*(rRNA)", "\\1", stats1$stitle)
    stats1$stitle <- gsub("(rRNA).*", "\\1", stats1$stitle)
    stats1$abundance2 <- stats1$abundance
    stats1$abundance[stats1$stitle == "ribosomal_RNA"] <- NA
    stats1$abundance[stats1$stitle == "rRNA"] <- NA
    stats1 <- subset(stats1, select=c(3,4,5))
    stats1$uniqseq <- 0
    stats1$uniqseq[stats1$abundance >= 0] <- 1
    stats1$abundance <- ifelse(is.na(stats1$abundance), stats1$abundance2, stats1$abundance)
    stats1 <- subset(stats1, select=-c(abundance2))
    stats2 <- ddply(stats1, c("phylum","uniqseq"), summarise, Nt = length(phylum))
    stats3 <- as.data.frame(table(stats2$phylum)); colnames(stats3)[1] <- c("phylum")
    stats1 <- merge(stats1, stats3, by=c("phylum"), all.x=TRUE)
    stats1$abundance[stats1$uniqseq == "0" & stats1$Freq == "2" ] <- NA
    stats1 <- ddply(stats1, c("phylum"), summarise,
                    N    = length(phylum),
                    mean = mean(abundance, na.rm=TRUE),
                    sd   = sd(abundance, na.rm=TRUE)
    )
    stats2 <- merge(stats2, stats3, by=c("phylum"), all.x=TRUE)
    stats2$Nt[stats2$uniqseq == "0" & stats2$Freq == "2" ] <- NA; stats2 <- na.omit(stats2)
    stats1 <- merge(stats1, stats2, by=c("phylum"), all.x=TRUE)
    stats1$se <- stats1$sd/stats1$Nt    
    stats1 <- subset(stats1, select=c(1,3,2,8))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

