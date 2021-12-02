#### PLOTS FOR MANUSCRIPT

library(ggplot2)

#######################################
############# FOLDED SFS ##############
#######################################

folded_sfs <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/MAF_Bins.txt",
                         sep="\t", header = T)
# order by severity
folded_sfs$Annotation <- factor(folded_sfs$Annotation, levels=c("NonCoding", "Synonymous",
                                                                "Tolerated", "AllDel", "StopLostGained"))

p <- ggplot(folded_sfs, aes(x=breaks, y=prop, fill=Annotation))
p + geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c("#999999", "#00AFBB", "#C3D7A4", "#FC4E07", "#E7B800")) +
  ylab("Proportion") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
                     labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", 
                                "[0.20,0.25)", "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)"))

#######################################
############ UNFOLDED SFS #############
#######################################

unfolded_sfs1 <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/DerivedFreq_Bins.txt",
                         sep="\t", header = T)

unfolded_sfs2 <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/Derived_dSNP_freqbins.txt",
                            sep="\t", header = T)

# use the dSNP annotation for deleterious
unfolded_sfs <- rbind.data.frame(unfolded_sfs1[which(unfolded_sfs1$Annotation!="AllDel"),], unfolded_sfs2)
# order by severity
unfolded_sfs$Annotation <- factor(unfolded_sfs$Annotation, levels=c("NonCoding", "Synonymous",
                                                                "Tolerated", "dSNP", "StopLostGained"))

p2 <- ggplot(unfolded_sfs, aes(x=breaks, y=prop, fill=Annotation))
p2 + geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c("#999999", "#00AFBB", "#C3D7A4", "#FC4E07", "#E7B800")) +
  ylab("Proportion") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                                0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00),
                     labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", "[0.20,0.25)", 
                                "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)",
                                "[0.50, 0.55)", "[0.55, 0.60)", "[0.60, 0.65)", "[0.65,0.70)", "[0.70, 0.75)", 
                                "[0.75, 0.80)", "[0.80, 0.85)", "[0.85, 0.90)", "[0.90, 0.95)", "[0.95, 1.00)"))


