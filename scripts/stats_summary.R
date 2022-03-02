#!/usr/bin/env Rscript


args <- commandArgs(TRUE)
sighits <- args[1]
minUniqRead <- args[2]
libdir <- args[4]

.libPaths( c( .libPaths(), libdir) )
library(plyr)
library(data.table)




if (args[3] == "strain"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,2,8))
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
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,staxids=stats1$staxids), FUN=sum)
    stats1 <- subset(stats1, select=c(3,2))
    names(stats1)[1] <- "abundance"
    # stats1 <- subset(stats1, stats1$abundance > 1)
    stats1 <- ddply(stats1, c("staxids"), summarise,
                    N    = length(staxids),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
    stats1 <- subset(stats1, stats1$uniq_reads >= minUniqRead)
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[3] == "species"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,3,13))
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
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,species=stats1$species), FUN=sum)
    stats1 <- subset(stats1, select=c(3,2))
    names(stats1)[1] <- "abundance"
    stats1 <- ddply(stats1, c("species"), summarise,
                    N    = length(species),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
    stats1 <- subset(stats1, stats1$uniq_reads >= minUniqRead)
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[3] == "genus"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,3,14))
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
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,genus=stats1$genus), FUN=sum)
    stats1 <- subset(stats1, select=c(3,2))
    names(stats1)[1] <- "abundance"
    stats1 <- ddply(stats1, c("genus"), summarise,
                    N    = length(genus),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
    stats1 <- subset(stats1, stats1$uniq_reads >= minUniqRead)
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[3] == "family"){
fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,3,15))
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
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,family=stats1$family), FUN=sum)
    stats1 <- subset(stats1, select=c(3,2))
    names(stats1)[1] <- "abundance"
    stats1 <- ddply(stats1, c("family"), summarise,
                    N    = length(family),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
    stats1 <- subset(stats1, stats1$uniq_reads >= minUniqRead)
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[3] == "order"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,3,16))
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
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,order=stats1$order), FUN=sum)
    stats1 <- subset(stats1, select=c(3,2))
    names(stats1)[1] <- "abundance"
    stats1 <- ddply(stats1, c("order"), summarise,
                    N    = length(order),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
    stats1 <- subset(stats1, stats1$uniq_reads >= minUniqRead)
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[3] == "class"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,3,17))
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
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,class=stats1$class), FUN=sum)
    stats1 <- subset(stats1, select=c(3,2))
    names(stats1)[1] <- "abundance"
    stats1 <- ddply(stats1, c("class"), summarise,
                    N    = length(class),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
    stats1 <- subset(stats1, stats1$uniq_reads >= minUniqRead)
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}

if (args[3] == "phylum"){
  fileName <- c(sighits)
  stats1 <- read.delim(file=fileName, header=T, sep="\t", fill= T, quote="", check.names = T)
  if (nrow(stats1) > 0) {
    stats1 <- subset(stats1, select=c(1,3,18))
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
    stats1 <- aggregate(stats1$abundance, by=list(sseqid=stats1$sseqid,phylum=stats1$phylum), FUN=sum)
    stats1 <- subset(stats1, select=c(3,2))
    names(stats1)[1] <- "abundance"
    stats1 <- ddply(stats1, c("phylum"), summarise,
                    N    = length(phylum),
                    mean = mean(abundance),
                    sd   = sd(abundance),
                    se   = sd / sqrt(N)
    )
    stats1 <- subset(stats1, select=c(1,3,2,5))
    colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")
    stats1$rel_stderr <- (stats1$stderr/stats1$mean)*100
    stats1 <- subset(stats1, stats1$uniq_reads >= minUniqRead)
  }
  if (nrow(stats1) == 0) { stats1 <- NULL; stats1 <- data.frame(matrix(nrow=0, ncol = 4)); colnames(stats1) <- c("taxid","mean","uniq_reads","stderr")}
  write.table(stats1,"stats1.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
}
