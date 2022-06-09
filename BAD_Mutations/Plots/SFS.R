## Site Frequency Spectrum Figures

source("BAD_Mutations/Variant_analyses/Functions.R")


#######################################
###### VARIANT CLASSES - FOLDED #######
#######################################

folded_sfs <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/MAF_Bins.txt",
                                sep="\t", header = T)

str(folded_sfs)

# order by severity
folded_sfs$Annotation <- factor(folded_sfs$Annotation, levels=c("NonCodingNodups", "SynonymousNodups",
                                                                "Tolerated", "MissenseOther", "AllDel", 
                                                                "StartStopLostGainednoDups", "SpliceAcceptorDonorNodups"))

AnnotationsToKeep <- c("SynonymousNodups", "Tolerated","AllDel")
folded_sfs_subset <- subset(folded_sfs, Annotation %in% AnnotationsToKeep)

p <- ggplot(folded_sfs_subset, aes(x=breaks, y=prop, fill=Annotation))
#p <- ggplot(folded_sfs, aes(x=breaks, y=prop, fill=Annotation))
p + geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c(Col_Synonymous, Col_Tolerated, Col_Deleterious),
                    name = "Variant Class", labels = c("Synonymous", "Tolerated", "Deleterious")) +
  ylab("Proportion") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
                     labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", 
                                "[0.20,0.25)", "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)")) +
  theme(panel.grid.minor = element_blank())

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Figure1_FoldedSFS.pdf",
       width = 10, height = 5)


#######################################
##### VARIANT CLASSES - UNFOLDED ######
#######################################


