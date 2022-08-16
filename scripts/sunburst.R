#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
perc <- args[2]
nlayers <- unlist(strsplit(args[3],","))
layers <- args[3]
libdir <- args[4]
taxlevel <- args[5]

.libPaths( c( .libPaths(), libdir) )
library(dplyr, quietly=T)
library(stringi, quietly=T)
library(plotme, quietly=T)
library(htmlwidgets, quietly=T)


perc <- as.numeric(perc)
sunburst <- read.delim(args[1], header=T, sep="\t", check.names=FALSE, quote="", fill=TRUE)
if (taxlevel == "strain") {
  sunburst <- subset(sunburst, select=-c(tax_id,kingdom,domain))
} else {
  sunburst <- subset(sunburst, select=-c(kingdom,domain))
}
colnames(sunburst) <- gsub("_mean","",colnames(sunburst))
sunburst <- data.frame(sunburst)
for (k in 1:(ncol(sunburst)-7)){
  sunburst[,k] <- as.numeric(as.character(sunburst[,k]))
}
sunburst$percent <- ((rowSums(sunburst[,1:(ncol(sunburst)-7)] > "0"))/(ncol(sunburst)-7))*100
sunburst <- subset(sunburst, percent >= perc)
sunburst <- subset(sunburst, select=-c(percent))
remtaxa <- nrow(sunburst)
# if (taxlevel == "strain") {
#   sunburst_virus <- sunburst[grepl("Viruses", sunburst$phylum),]
#   sunburst <- sunburst[!grepl("viricota", sunburst$phylum),]
#   sunburst_virus$species <- sunburst_virus$taxname
#   sunburst_virus_genus <- sunburst_virus[is.na(sunburst_virus$genus),]
#   sunburst_virus <- sunburst_virus[!(is.na(sunburst_virus$genus)),]
#   sunburst_virus_genus$genus <- sub("\\_.*","",sunburst_virus_genus$taxname)
#   sunburst_virus <- rbind(sunburst_virus, sunburst_virus_genus)
#   sunburst$genus <- sub("\\_.*","",sunburst$taxname)
#   sunburst$species <- sub("\\_.*","",(sub(".*?_","",sunburst$taxname)))
#   sunburst <- rbind(sunburst, sunburst_virus)
# }

sunburst[["phylum"]][is.na(sunburst[["kingdom"]])] <- "Phylum Unclassified"
sunburst[["class"]][is.na(sunburst[["kingdom"]])] <- "Class Unclassified"
sunburst[["order"]][is.na(sunburst[["kingdom"]])] <- "Order Unclassified"
sunburst[["family"]][is.na(sunburst[["kingdom"]])] <- "Family Unclassified"
sunburst[["genus"]][is.na(sunburst[["kingdom"]])] <- "Genus Unclassified"
# if (taxlevel == "strain") {
#   for (i in 1:nrow(sunburst)){
#     if(is.na(sunburst$species[i])){
#       sunburst$species[i] <- sunburst$taxname[i]
#     }
#   }
# }
# if (taxlevel == "species") {
#   for (i in 1:nrow(sunburst)){
#     if(is.na(sunburst$species[i])){
#       sunburst$species[i] <- sunburst$taxname[i]
#     }
#   }
# }
sunburst$average <- (rowSums(sunburst[,1:(ncol(sunburst)-7)])/(ncol(sunburst)-7))

layers <- gsub(",","_",layers)
dir.create("sunburst", showWarnings = FALSE)
rel_abun <- sunburst %>%
  count(sunburst[,(match(nlayers,names(sunburst)))], wt=average) %>%
  count_to_sunburst()
htmlwidgets::saveWidget(rel_abun, paste("./sunburst/taxa_profile_rel_abundance_",perc,"perc_",remtaxa,"_taxa_",layers,".html",sep=""), selfcontained=F)

diversity <- sunburst %>%
  count(sunburst[,(match(nlayers,names(sunburst)))],wt=NULL) %>%
  count_to_sunburst()
htmlwidgets::saveWidget(diversity, paste("./sunburst/taxa_profile_diversity_",perc,"perc_",remtaxa,"_taxa_",layers,".html",sep=""), selfcontained=F)


