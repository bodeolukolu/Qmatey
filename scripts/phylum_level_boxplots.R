#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
dfm <- read.delim(args[1], header=T, sep="\t", check.names =FALSE, fill=TRUE)
dfrm <- read.delim(args[2], header=T, sep="\t", check.names =FALSE, fill=TRUE)
dfu <- read.delim(args[3], header=T, sep="\t", check.names = FALSE, fill=TRUE)
dfe <- read.delim(args[4], header=T, sep="\t", check.names = FALSE, fill= TRUE)
dfre <- read.delim(args[5], header=T, sep="\t", check.names = FALSE, fill= TRUE)
perc <- as.numeric(args[6])
libdir <- args[7]

dfm$percent <- (rowSums(dfm[,2:(ncol(dfm)-2)] > "0")/(ncol(dfm)-3))*100 
dfm <- subset(dfm, percent >= perc)
dfm <- subset(dfm, select=-c(percent))
remtaxa <- nrow(dfm)
dfrm <- dfrm[c(rownames(dfm)),]
dfu <- dfu[c(rownames(dfm)),]
dfe <- dfe[c(rownames(dfm)),]
dfre <- dfre[c(rownames(dfm)),]


.libPaths( c( .libPaths(), libdir) )
library(dplyr, quietly=T)
library(plotly, quietly=T)
library(htmlwidgets, quietly=T)



data.pipe<-function(dfm,dfrm,dfu,dfe,dfre,cut){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  vDrop<-c("kingdom","domain")
  mirror<-dfm[,(which(colnames(dfm) %!in% vDrop))]
  mirror$phylum = as.character(mirror$phylum)
  mirror <- mirror[c(2:ncol(mirror),1,1)]
  colnames(mirror)[ncol(mirror)] <- c("Phylum")
  
  #set up empty arrays for the data values
  rowNum<-0
  taxid<-c()
  means<-c()
  uniq<-c()
  err<-c()
  rerr<-c()
  phylum<-c()
  lines<-0
  #loop to input data into selected arrays
  for (row in 1:nrow(mirror) ) {
    lines=lines+1
    for (col in 1:(ncol(mirror)-2)) {
      rowNum<- rowNum + 1
      taxid[rowNum]<-mirror[row,ncol(mirror)-1]
      phylum[rowNum]<-mirror[row,ncol(mirror)]
      means[rowNum]<-mirror[row,col]
      rmeans[rowNum]<-mirror[row,col]
      uniq[rowNum]<-dfu[row,col+1]
      err[rowNum]<-dfe[row,col+1]
      rerr[rowNum]<-dfre[row,col+1]
    }
    
  }
  
  #create dataframe to hold the data
  reform<-data.frame(matrix(nrow=rowNum, ncol = 6))
  
  #convert the arrays to vectors 
  taxid<-as.vector(unlist(taxid))
  means<-as.vector(unlist(means))
  rmeans<-as.vector(unlist(rmeans))
  uniq<-as.vector(unlist(uniq))
  err<-as.vector(unlist(err))
  rerr<-as.vector(unlist(rerr))
  phy<-as.vector(unlist(phylum))
  #set variables equal to corresponding vectors of data
  reform$X1<-as.factor(taxid)
  reform$X2<-means
  reform$X3<-rmeans
  reform$X4<-uniq
  reform$X5<-err
  reform$X6<-rerr
  reform$X7<-phy
  reform<-subset(reform, reform$X2 >0)
  reform<-subset(reform, reform$X3 >0)
  reform<-subset(reform, reform$X4 >0)
  reform<-subset(reform, reform$X5 >0)
  reform<-subset(reform, reform$X6 >0)
  #rename the variables
  colnames(reform)<-c("taxid","covMean","Unique Reads","errors","rel_errors","Phylum")
  colnames(reform)<-c("taxid","covMean","RelMean","Unique Reads","errors","rel_errors","Phylum")
  #create total coverage data by multiplying the coverage by the unique reads
  reform$'Unique Reads' <- as.numeric(as.character(reform$'Unique Reads'))
  reform$'covMean' <- as.numeric(as.character(reform$'covMean'))
  reform$'RelMean' <- as.numeric(as.character(reform$'RelMean'))
  reform$'errors' <- as.numeric(as.character(reform$'errors'))
  reform$'rel_errors' <- as.numeric(as.character(reform$'rel_errors'))
  reform$'Tcov'<-reform$'covMean'*reform$'Unique Reads'
  
  #cut data that makes up small percent of data
  met<-dfm
  met<-subset.data.frame(met, select=-c(kingdom,domain))
  met <- met[c(2:(ncol(met)-1),1,ncol(met))]
  met<-subset.data.frame(met, select=c(ncol(met),1:(ncol(met)-1)))
  met$count <-rowSums(met[,2:ncol(met)] != "0")
  met$percent <-((met$count)/(ncol(met)-2))*100
  bxD<-subset.data.frame(met, met$percent>= (5)/100)
  reform<-subset(reform, reform$taxid %in% bxD$phylum)
  for (r in c(2,3,4,6)) { 
    reform[,r][reform[,r]==0]<-NA
  }
  
  
  return(reform)
}

reform <- data.pipe(dfm, dfrm, dfu, dfe, dfre, cut)
df <- reform


uniq_box<-function(df, dfu){
  uniq <- subset(dfu, select=-c(kingdom,domain))
  uniq <- uniq[c(2:(ncol(uniq)-1),1,ncol(uniq))]
  fence<-data.frame(matrix(nrow=nrow(uniq), ncol = 2))
  taxa<-c()
  up<-c()
  for (pl in 1:nrow(uniq)) {
    taxa[pl]<-uniq[pl,ncol(uniq)]
    vec<-as.numeric(uniq[pl,2:(ncol(uniq)-1)])
    iqr<-IQR(vec, na.rm = TRUE)
    # UpperInner Fence
    up[pl] <- quantile(vec, .75, na.rm = TRUE) + (1.5*iqr)
  }
  #set variables equal to corresponding vectors of data
  fence$X1<-taxa
  fence$X2<-up
  #remove
  colnames(fence)<-c("taxid","Upper")
  #make the boxplot
  xlim<-max(fence$Upper)
  reform$Phylum[is.na(reform$Phylum)] = "Unknown"
  u<-plot_ly(reform, x = ~`Unique Reads`, y = ~taxid, color = ~Phylum, type = "box", boxpoints = FALSE)%>%layout(title = paste("Boxplot of Unique Reads (",remtaxa," taxa)",sep=""), xaxis = list(range=c(0,xlim), title = "Unique Reads"), yaxis = list(size = 1, title = "Taxaname"))
  htmlwidgets::saveWidget(as_widget(u), paste("phylum_level_unique_reads_",perc,"perc.html",sep=""), selfcontained=FALSE)
}
uniq_box(df, dfu)

mean_box<-function(df, dfm){
  mean <- subset(dfm, select=-c(kingdom,domain))
  mean <- mean[c(2:(ncol(mean)-1),1,ncol(mean))]
  fence<-data.frame(matrix(nrow=nrow(mean), ncol = 2))
  taxa<-c()
  up<-c()
  for (pl in 1:nrow(mean)) {
    taxa[pl]<-mean[pl,ncol(mean)]
    vec<-as.numeric(mean[pl,2:(ncol(mean)-1)])
    iqr<-IQR(vec, na.rm = TRUE)
    # UpperInner Fence
    up[pl] <- quantile(vec, .75, na.rm = TRUE) + (1.5*iqr)
  }
  #set variables equal to corresponding vectors of data
  fence$X1<-taxa
  fence$X2<-up
  #remove
  colnames(fence)<-c("taxid","Upper")
  #make the boxplot
  xlim<-max(fence$Upper)
  reform$Phylum[is.na(reform$Phylum)] = "Unknown"
  u<-plot_ly(reform, x = ~`covMean`, y = ~taxid, color = ~Phylum, type = "box", boxpoints = FALSE)%>%layout(title = paste("Boxplot of Mean Absolute Abundance (",remtaxa," taxa)",sep=""), xaxis = list(range=c(0,xlim), title = "Mean Reads"), yaxis = list(size = 1, title = "Taxaname"))
  htmlwidgets::saveWidget(as_widget(u), paste("phylum_level_mean_reads_",perc,"perc.html",sep=""), selfcontained=FALSE)
}
mean_box(df, dfm)

mean_box<-function(df, dfrm){
  mean <- subset(dfrm, select=-c(kingdom,domain))
  mean <- mean[c(1:(ncol(mean)-1),ncol(mean),ncol(mean))]
  fence<-data.frame(matrix(nrow=nrow(mean), ncol = 2))
  taxa<-c()
  up<-c()
  for (pl in 1:nrow(mean)) {
    taxa[pl]<-mean[pl,ncol(mean)]
    vec<-as.numeric(mean[pl,2:(ncol(mean)-1)])
    iqr<-IQR(vec, na.rm = TRUE)
    # UpperInner Fence
    up[pl] <- quantile(vec, .75, na.rm = TRUE) + (1.5*iqr)
  }
  #set variables equal to corresponding vectors of data
  fence$X1<-taxa
  fence$X2<-up
  #remove
  colnames(fence)<-c("taxid","Upper")
  #make the boxplot
  xlim<-max(fence$Upper)
  reform$Phylum[is.na(reform$Phylum)] = "Unknown"
  u<-plot_ly(reform, x = ~`covMean`, y = ~taxid, color = ~Phylum, type = "box", boxpoints = FALSE)%>%layout(title = paste("Boxplot of Mean Relative Abundance (",remtaxa," taxa)",sep=""), xaxis = list(range=c(0,xlim), title = "Mean Reads"), yaxis = list(size = 1, title = "Taxaname"))
  htmlwidgets::saveWidget(as_widget(u), paste("phylum_level_rel_mean_reads_",perc,"perc.html",sep=""), selfcontained=FALSE)
}
mean_box(df, dfrm)

error_box<-function(df, dfe){
  error <- subset(dfe, select=-c(kingdom,domain))
  error <- error[c(2:(ncol(error)-1),1,ncol(error))]
  fence<-data.frame(matrix(nrow=nrow(error), ncol = 2))
  taxa<-c()
  up<-c()
  for (pl in 1:nrow(error)) {
    taxa[pl]<-error[pl,ncol(error)]
    vec<-as.numeric(error[pl,2:(ncol(error)-1)])
    iqr<-IQR(vec, na.rm = TRUE)
    # UpperInner Fence
    up[pl] <- quantile(vec, .75, na.rm = TRUE) + (1.5*iqr)
  }
  #set variables equal to corresponding vectors of data
  fence$X1<-taxa
  fence$X2<-up
  #remove
  colnames(fence)<-c("taxid","Upper")
  #make the boxplot
  xlim<-max(fence$Upper)
  reform$Phylum[is.na(reform$Phylum)] = "Unknown"
  u<-plot_ly(reform, x = ~errors, y = ~taxid, color = ~Phylum, type = "box", boxpoints = FALSE)%>%layout(title= paste("Boxplot of Standard Error (",remtaxa," taxa)",sep=""), xaxis = list(titlerange=c(0,xlim), title = "Standard Error"), yaxis = list(size = 1, title = "Taxaname"))
  htmlwidgets::saveWidget(as_widget(u), paste("phylum_level_errors_",perc,"perc.html",sep=""), selfcontained=FALSE)
}
error_box(df, dfe)

rel_error_box<-function(df, dfre){
  rel_error <- subset(dfre, select=-c(kingdom,domain))
  rel_error <- rel_error[c(2:(ncol(rel_error)-1),1,ncol(rel_error))]
  fence<-data.frame(matrix(nrow=nrow(rel_error), ncol = 2))
  taxa<-c()
  up<-c()
  for (pl in 1:nrow(rel_error)) {
    taxa[pl]<-rel_error[pl,ncol(rel_error)]
    vec<-as.numeric(rel_error[pl,2:(ncol(rel_error)-1)])
    iqr<-IQR(vec, na.rm = TRUE)
    # UpperInner Fence
    up[pl] <- quantile(vec, .75, na.rm = TRUE) + (1.5*iqr)
  }
  #set variables equal to corresponding vectors of data
  fence$X1<-taxa
  fence$X2<-up
  #remove
  colnames(fence)<-c("taxid","Upper")
  #make the boxplot
  xlim<-max(fence$Upper)
  reform$Phylum[is.na(reform$Phylum)] = "Unknown"
  u<-plot_ly(reform, x = ~errors, y = ~taxid, color = ~Phylum, type = "box", boxpoints = FALSE)%>%layout(title= paste("Boxplot of Relative Standard Error (",remtaxa," taxa)",sep=""), xaxis = list(titlerange=c(0,xlim), title = "Relative Standard Error (% RSE)"), yaxis = list(size = 1, title = "Taxaname"))
  htmlwidgets::saveWidget(as_widget(u), paste("phylum_level_rel_errors_",perc,"perc.html",sep=""), selfcontained=FALSE)
}
rel_error_box(df, dfre)

