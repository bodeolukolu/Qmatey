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
if (taxlevel == "strain") { sunburst$percent <- ((rowSums(sunburst[,1:(ncol(sunburst)-7)] > "0"))/(ncol(sunburst)-7))*100 }
if (taxlevel == "species") { sunburst$percent <- ((rowSums(sunburst[,1:(ncol(sunburst)-6)] > "0"))/(ncol(sunburst)-6))*100 }
if (taxlevel == "genus") { sunburst$percent <- ((rowSums(sunburst[,1:(ncol(sunburst)-5)] > "0"))/(ncol(sunburst)-5))*100 }
if (taxlevel == "family") { sunburst$percent <- ((rowSums(sunburst[,1:(ncol(sunburst)-4)] > "0"))/(ncol(sunburst)-4))*100 }
if (taxlevel == "order") { sunburst$percent <- ((rowSums(sunburst[,1:(ncol(sunburst)-3)] > "0"))/(ncol(sunburst)-3))*100 }
if (taxlevel == "class") { sunburst$percent <- ((rowSums(sunburst[,1:(ncol(sunburst)-2)] > "0"))/(ncol(sunburst)-2))*100 }
if (taxlevel == "phylum") { sunburst$percent <- ((rowSums(sunburst[,1:(ncol(sunburst)-1)] > "0"))/(ncol(sunburst)-1))*100 }
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

if (taxlevel == "genus") {
  sunburst[["phylum"]][is.na(sunburst[["kingdom"]])] <- "Phylum Unclassified"
  sunburst[["class"]][is.na(sunburst[["kingdom"]])] <- "Class Unclassified"
  sunburst[["order"]][is.na(sunburst[["kingdom"]])] <- "Order Unclassified"
  sunburst[["family"]][is.na(sunburst[["kingdom"]])] <- "Family Unclassified"
  sunburst[["genus"]][is.na(sunburst[["kingdom"]])] <- "Genus Unclassified"
}
if (taxlevel == "family") {
  sunburst[["phylum"]][is.na(sunburst[["kingdom"]])] <- "Phylum Unclassified"
  sunburst[["class"]][is.na(sunburst[["kingdom"]])] <- "Class Unclassified"
  sunburst[["order"]][is.na(sunburst[["kingdom"]])] <- "Order Unclassified"
  sunburst[["family"]][is.na(sunburst[["kingdom"]])] <- "Family Unclassified"
}
if (taxlevel == "order") {
  sunburst[["phylum"]][is.na(sunburst[["kingdom"]])] <- "Phylum Unclassified"
  sunburst[["class"]][is.na(sunburst[["kingdom"]])] <- "Class Unclassified"
  sunburst[["order"]][is.na(sunburst[["kingdom"]])] <- "Order Unclassified"
}
if (taxlevel == "class") {
  sunburst[["phylum"]][is.na(sunburst[["kingdom"]])] <- "Phylum Unclassified"
  sunburst[["class"]][is.na(sunburst[["kingdom"]])] <- "Class Unclassified"
}
if (taxlevel == "phylum") {
  sunburst[["phylum"]][is.na(sunburst[["kingdom"]])] <- "Phylum Unclassified"
}

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
if (taxlevel == "strain") { sunburst$average <- (rowSums(sunburst[,1:(ncol(sunburst)-7)])/(ncol(sunburst)-7)) }
if (taxlevel == "species") { sunburst$average <- (rowSums(sunburst[,1:(ncol(sunburst)-6)])/(ncol(sunburst)-6)) }
if (taxlevel == "genus") { sunburst$average <- (rowSums(sunburst[,1:(ncol(sunburst)-5)])/(ncol(sunburst)-5)) }
if (taxlevel == "family") { sunburst$average <- (rowSums(sunburst[,1:(ncol(sunburst)-4)])/(ncol(sunburst)-4)) }
if (taxlevel == "order") { sunburst$average <- (rowSums(sunburst[,1:(ncol(sunburst)-3)])/(ncol(sunburst)-3)) }
if (taxlevel == "class") { sunburst$average <- (rowSums(sunburst[,1:(ncol(sunburst)-2)])/(ncol(sunburst)-2)) }
if (taxlevel == "phylum") { sunburst$average <- (rowSums(sunburst[,1:(ncol(sunburst)-1)])/(ncol(sunburst)-1)) }

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


