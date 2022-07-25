### PCA plots

source("BAD_Mutations/Variant_analyses/Functions.R")

library(ggplot2)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(scales)

###############################
###### SET UP DATAFRAME #######
###############################

smartpca <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/Sunflower_SAM.pca.evec",
                       col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5",
                                   "PC6", "PC7", "PC8", "PC9", "PC10", "Pop"))
# 30 outliers removed

# no outliers removed:
#smartpca <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/Sunflower_SAM_all.pca.evec",
#                       col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5",
#                                   "PC6", "PC7", "PC8", "PC9", "PC10", "Pop"))
# extreme outliers: PPN285(Hopi), PPN083 (landrace), PPN093 (landrace)
# PPN204, PPN288, PPN203, PPN217 (all RHA-Oil)
# PPN150, PPN019, PPN090 (HA-Oil)
# PPN170 (HA-NonOiil), PPN079 (introgressed)

# combine with sample info
key <- read.csv("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Sunflower_Mutation_Load/BAD_Mutations/LineKeywINFO.csv", 
                header=T)

pca_wInfo <- merge(smartpca, key, by.x="Sample", by.y = "PPN")

smartpca[which(!smartpca$Sample %in% pca_wInfo$Sample),] # PI 531071 and SF-33

levels(pca_wInfo$Groups)
pca_wInfo$Groups <- factor(pca_wInfo$Groups,
                              levels=c("HA-Oil", "RHA-Oil",
                                       "HA-NonOil", "RHA-NonOil",
                                       "introgressed", "landrace_OPV"))

# make groups for graphing
#pca_wInfo$group <- as.factor(ifelse(pca_wInfo$heterotic_group=="OPV" |
#                            pca_wInfo$heterotic_group=="landrace",
#                          as.character(pca_wInfo$heterotic_group), 
#                          as.character(pca_wInfo$Groups)))

###############################
######### MAKE PLOT ###########
###############################

show_col(pal_jco("default")(10))

HA_Oil <- c("#EFC000FF")
RHA_Oil <- c("#8F7700FF")
HA_NonOil <- c("#7AA6DCFF")
RHA_NonOil <- c("#4A6990FF")
#RHA_NonOil <- c("#003C67FF")
#RHA_NonOil <- c("#0073C2FF")
introgressed <- c("#CD534CFF")
landrace_OPV <- c("#A73030FF")

ggplot(data = pca_wInfo, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=Groups), size=3) +
  scale_color_jco() +
  theme_minimal()

ggplot(data = pca_wInfo, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill=Groups, shape=Groups, color=Groups), size=4, alpha=0.8) +
  scale_shape_manual(values=c(21, 21, 21, 21, 23, 8)) +
  scale_fill_manual(values=c(HA_Oil, RHA_Oil, HA_NonOil, RHA_NonOil,
                              introgressed, landrace_OPV)) +
  scale_color_manual(values=c("black", "black", "black", "black", "black", landrace_OPV)) +
  #geom_text(aes(label=SequenceName)) +
  theme_minimal()

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/PCA.pdf",
       width = 10, height = 8)

###############################
##### PLOT FOR CHROM 10 #######
###############################

pca10 <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/Sunflower_SAM10.pca.evec",
                       col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5",
                                   "PC6", "PC7", "PC8", "PC9", "PC10", "Pop"))

pca10_wInfo <- merge(pca10, key, by.x="Sample", by.y = "PPN")

# subset to only HA and RHA
pca10_hetgroups <- subset(pca10_wInfo, Groups!="landrace_OPV" &
                            Groups!="introgressed")

ggplot(data = pca10_hetgroups, aes(x=PC1, y=PC2)) +
ggplot(data = pca10_wInfo, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=Groups), size=3, alpha=0.8) +
  scale_color_jco() +
  theme_minimal() +
  geom_text(aes(label=SequenceName, color=Groups))

