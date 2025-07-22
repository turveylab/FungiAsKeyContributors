library(tidyverse)
library(phyloseq)
library(ggpubr)
library(ggforce)
library(patchwork)
library(phytools)
library(ape)

#load taxa table
taxa.table<-read_tsv("~/Documents/MetagenomicsTurvey/Sequence_outputs/ITS2_Qiime/taxonomy.tsv")
taxa.table<-subset(taxa.table, !is.na(Phylum))
taxa.table<-subset(taxa.table, !is.na(Genus))
taxa.table<-subset(taxa.table,Genus != 'NA unknown' )
taxa.table$Class<-if_else(is.na(taxa.table$Class), paste(taxa.table$Phylum,"unknown"),taxa.table$Class)
taxa.table$Order<-if_else(is.na(taxa.table$Order), paste(taxa.table$Class,"unknown"),taxa.table$Order)
taxa.table$Order <-str_replace(taxa.table$Order,pattern = "unknown unknown", replace="unknown")
taxa.table$Family<-if_else(is.na(taxa.table$Family), paste(taxa.table$Order,"unknown"),taxa.table$Family)
taxa.table$Family <-str_replace(taxa.table$Family,pattern = "unknown unknown", replace="unknown")
taxa.table$Genus<-if_else(is.na(taxa.table$Genus), paste(taxa.table$Family,"unknown"),taxa.table$Genus)
taxa.table$Genus <-str_replace(taxa.table$Genus,pattern = "unknown unknown", replace="unknown")
taxa.table$Species<-if_else(is.na(taxa.table$Species), paste(taxa.table$Genus,"unknown"),taxa.table$Species)
taxa.table$Species <-str_replace(taxa.table$Species,pattern = "unknown unknown", replace="unknown")

otu.table<-otu.table[rownames(otu.table) %in% taxa.table$`Feature ID`,]

metadata<-read_csv("~/Documents/MetagenomicsTurvey/MetadataFiles/CHILD_metagenomics_metadata_all.csv")

write_tsv(otu.table,"~/Documents/MetagenomicsTurvey/Sequence_outputs/ITS2_Qiime/its2-nounk-otu-table.tsv")
write_tsv(taxa.table,"~/Documents/MetagenomicsTurvey/Sequence_outputs/ITS2_Qiime/its2-nounk-taxa-table.tsv")
write.tree(pruned.tree, "~/Documents/MetagenomicsTurvey/Sequence_outputs/ITS2_Qiime/its2-nounk.nwk")
#xlim = c(100000,5000000)
#sequencing length compare

otu.table.pr.4k<-colSums(otu.table) >= 4000
otu.table.pr.4k<-otu.table[,otu.table.pr.4k]

sequencedepth_filt.pr_4k <- colSums(otu.table.pr.4k[sapply(otu.table.pr.4k, is.numeric)])
hist(sequencedepth_filt.pr_4k, breaks=1000, col="grey", main="", ylab="Number of Samples", xlab="Sequencing Depth", xlim = c(1000,50000))

colnames(otu.table.pr.4k) <- gsub('\\_', '.', colnames(otu.table.pr.4k))
taxa.table <- as.data.frame(taxa.table)
rownames(taxa.table) <- taxa.table$`Feature ID`
taxa.table <- taxa.table[ , -1]
metadata <- as.data.frame(metadata)
metadata$Barcodes <- gsub('\\-', '.', metadata$Barcodes)
rownames(metadata) <- metadata$Barcodes
otu.table.pr.4k.t <- t(otu.table.pr.4k)
metadata <- subset(metadata, rownames(metadata)%in%rownames(otu.table.pr.4k.t))
reorder_idx <- match(rownames(otu.table.pr.4k.t),rownames(metadata))
metadata <- metadata[reorder_idx, ]
otu.table.pr.4k.t<- subset(otu.table.pr.4k.t, rownames(otu.table.pr.4k.t)%in%rownames(metadata))
otu.table.pr.4k.t <- as.data.frame(t(otu.table.pr.4k.t))

OTU = otu_table(as.matrix(otu.table.pr.4k.t), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxa.table))
samples = sample_data(metadata)


ITS2_noUnk_4k <- phyloseq(OTU, TAX, samples)
ITS2_noUnk_4k

ITS2_noUnk_4k_genus <- tax_glom(ITS2_noUnk_4k, "Genus")
ITS2_noUnk_4k_genus
