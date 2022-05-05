#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
perc <- args[2]
mincorr <- unlist(strsplit(args[3],","))
maxcorr <- unlist(strsplit(args[4],","))
libdir <- args[5]

.libPaths( c( .libPaths(), libdir) )
library(ggplot2)
library(ggcorrplot)
library(reshape2)
library(gtools)
source(paste(libdir,"/CCLasso/R/cclasso.R",sep=""))


perc <- as.numeric(perc)
metag <- read.delim(args[1], header=T, sep="\t", check.names=FALSE, quote="", fill=TRUE)
metag <- subset(metag, select=-c(kingdom,domain))
metag <- subset(metag, select=c(ncol(metag),1:(ncol(metag)-1)))
metag <- setNames(data.frame(t(metag[,-1])), metag[,1])
metag <- data.frame(t(metag))
for (k in 1:ncol(metag)){
  metag[,k] <- as.numeric(as.character(metag[,k]))
}
metag$percent <- (rowSums(metag[,1:ncol(metag)] > "0")/ncol(metag))*100
metag <- subset(metag, percent >= perc)
metag <- subset(metag, select=-c(percent))
metag <- data.frame(t(metag))
taxa <- colnames(metag)

res_ccl_metag <- cclasso(x = metag, counts = T, k_cv=5, n_boot=100)
metag_corr <- res_ccl_metag$cor_w; row.names(metag_corr) <- taxa; colnames(metag_corr) <- taxa
metag_pmat <- res_ccl_metag$p_vals; row.names(metag_pmat) <- taxa; colnames(metag_pmat) <- taxa
axis_density <- 2000/nrow(metag_corr)
  

metag_corr[is.na(metag_corr)] <- 0
metag_pmat[is.na(metag_pmat)] <- 1
ordering <- hclust(dist(metag_corr), method = "average", members = NULL)
ordered <- as.data.frame(ordering[["labels"]]); ordered <- ordered[ordering[["order"]],]
metag_corr <- metag_corr[,ordered]; metag_corr <- metag_corr[ordered,]
metag_corr1 <- data.frame()
for ( j in 1:ncol(metag_corr)) {
  metag_corr0 <- as.data.frame(metag_corr[,j]); colnames(metag_corr0)[1] <- c("coeff")
  metag_corr0$Var1 <- colnames(metag_corr)[j]; metag_corr0$Var2 <- rownames(metag_corr0)
  metag_corr0 <- subset(metag_corr0, select=c(2,3,1))
  metag_corr1 <- rbind(metag_corr1, metag_corr0); metag_corr0 <- NULL
}
metag_pmat <- metag_pmat[,ordered]; metag_pmat <- metag_pmat[ordered,]
metag_pmat1 <- data.frame()
for ( j in 1:ncol(metag_pmat)) {
  metag_pmat0 <- as.data.frame(metag_pmat[,j]); colnames(metag_pmat0)[1] <- c("pvalue")
  metag_pmat0$Var1 <- colnames(metag_pmat)[j]; metag_pmat0$Var2 <- rownames(metag_pmat0)
  metag_pmat0 <- subset(metag_pmat0, select=c(2,3,1))
  metag_pmat1 <- rbind(metag_pmat1, metag_pmat0); metag_pmat0 <- NULL
}
corr <- cbind(metag_corr1, metag_pmat1)
corr <- subset(corr, select=-c(4:5))
write.table(corr,paste("Metagenome_",perc,"perc_phylum_corr.txt",sep=""), sep="\t",row.names=FALSE, col.names=TRUE, quote = F)

for (minc in (as.numeric(mincorr))) {
  for (maxc in (as.numeric(maxcorr))) {
    corr[,3][corr[,4] > 0.05] <- NA; corr[,4][corr[,4] > 0.05] <- NA
    corr[,4][corr[,3] < minc & corr[,3] > -maxc ] <- NA; corr[,3][corr[,3] < minc & corr[,3] > -maxc ] <- NA
    corr[,3][corr[,1] == corr[,2] ] <- NA
    corr <- na.omit(corr)
    corr$Var1 <- factor(corr$Var1, levels = ordered)
    corr$Var2 <- factor(corr$Var2, levels = ordered)
    ntaxa <- length(unique(c(corr[,1],corr[,2])))
    
    plot <- ggplot(corr, aes(x=Var1, y=Var2, fill=coeff, size= circle_size)) +
      geom_point(aes(size=-pvalue), shape=21) + scale_fill_gradient2(low="red", mid="white", high="cornflowerblue") +
      theme_bw() + coord_equal() + scale_size(range=c(1,10),guide = 'none') +
      labs(x="",y="",fill="Correlation\nCoefficient",size="p-value") +
      theme(axis.text.x=element_text(size=axis_density, angle=45, vjust=1, hjust=1, margin=margin(0,0,0,0)),
            axis.text.y=element_text(size=axis_density, margin=margin(0,0,0,0)), panel.grid.major=element_line(colour = "grey95"),
            legend.title=element_text(size=0.5*axis_density), legend.text=element_text(size=20),legend.key.size = unit(0.5, "in"),
            plot.title = element_text(size=axis_density)) +
      labs(title= paste("Phylum-Level Correlogram (",perc,"% missingness threshold, ",ntaxa," taxa)",sep=""))
    ggsave(filename=paste("Metagenome_phylum_",perc,"perc_","neg",maxc,"_pos",minc,"_corr.tiff",sep=""), plot=plot, width=35, height=35, dpi=600, compression = "lzw", limitsize = FALSE)
  }
}

