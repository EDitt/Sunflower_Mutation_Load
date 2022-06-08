#### PLOTS FOR MANUSCRIPT

library(ggplot2)
library(ggpubr)
library(viridis)
library(RColorBrewer)

#######################################
############# FOLDED SFS ##############
#######################################

folded_sfs <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/MAF_Bins.txt",
                         sep="\t", header = T)
#folded_sfs <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/MAF_Bins_noOneOrTwo.txt",
#                         sep="\t", header = T)

# order by severity
folded_sfs$Annotation <- factor(folded_sfs$Annotation, levels=c("NonCoding", "Synonymous",
                                                                "Tolerated", "AllDel", "StopLostGained"))

p <- ggplot(folded_sfs, aes(x=breaks, y=prop, fill=Annotation))
p + geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c("#999999", "#00AFBB", "#C3D7A4", "#FC4E07", "#E7B800"),
                    name = "Variant Class", labels = c("Non-coding", "Synonymous", "Tolerated", "Deleterious", "Stop Codon Lost/Gained")) +
  ylab("Proportion") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
                     labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", 
                                "[0.20,0.25)", "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)"))

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Figure1_FoldedSFS.pdf")

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
  scale_fill_manual(values= c("#999999", "#00AFBB", "#C3D7A4", "#FC4E07", "#E7B800"),
                    name = "Variant Class", labels = c("Non-coding", "Synonymous", "Tolerated", "Deleterious", "Stop Codon Lost/Gained")) +
  ylab("Proportion") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                                0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00),
                     labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", "[0.20,0.25)", 
                                "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)",
                                "[0.50, 0.55)", "[0.55, 0.60)", "[0.60, 0.65)", "[0.65,0.70)", "[0.70, 0.75)", 
                                "[0.75, 0.80)", "[0.80, 0.85)", "[0.85, 0.90)", "[0.90, 0.95)", "[0.95, 1.00)"))

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/FigureS1_UnfoldedSFS.pdf",
       width = 20, height = 5, units = "in")

# subset (only synonymous, tolerated, and dSNP)
unfolded_sfs_sub <- subset(unfolded_sfs, Annotation=="Synonymous" | Annotation=="Tolerated" | Annotation == "dSNP")

p3 <- ggplot(unfolded_sfs_sub, aes(x=breaks, y=prop, fill=Annotation))
p3 + geom_bar(stat="identity", position = position_dodge(), color="black") + 
  theme_minimal() +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_fill_viridis(discrete = TRUE, option = "D", alpha = 0.8)
  scale_fill_manual(values= c("#7570B3", "#1B9E77", "#D95F02"),
                    name = "Variant Class", labels = c("Synonymous", "Tolerated", "Deleterious")) +
  ylab("Proportion") +
  xlab("Frequency") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16, face = "bold"),
        panel.grid.minor = element_blank(),
        #legend.title = element_blank(),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = c(0.7, 0.5),
        legend.background = element_rect(fill = "white", size = 0.1, linetype = "solid")) +
scale_x_continuous(breaks = c(0.05, 0.25, 0.50, 0.75, 1.00),
                   labels = c("[0 - 0.05)", 
                              "[0.20 - 0.25)",
                              "[0.45 - 0.50)", 
                              "[0.70 - 0.75)", "[0.95 - 1.00)"))

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/UnfoldedSFS_subset.pdf",
       width = 20, height = 10, units = "in")

#######################################
##### RECOMBINATION & dSNP/CODON ######
#######################################

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/RecombinationDFs.RData")
# objects: Recomb_ChromCalcs (cM_Mbp info for all markers)

VariantNums <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/Derived_VariantNums.txt",
                          header=T, sep="\t")

VariantNums$MiddlePos <- VariantNums$StartPos + 5000000

# split by chromosome
VariantNums_Chromosome <- split(VariantNums, VariantNums$Chromosome)

colors <- c("dSNP per Codon" = "#FC4E07", "sSNP per Codon" = "#00AFBB",
            "cM/Mbp" = "blue")

Recomb_VarNum_Codon_plot <- function(variant_dataset, recombination_dataset, maxY, scaleNum) {
  plot <- ggplot(variant_dataset, aes(x=MiddlePos/1000000, y=(Number_dSNP/Num_codons))) +
    #geom_point(size=2, shape=24, fill=alpha("#FC4E07", 0.9)) +
    geom_smooth(method="loess", se=FALSE, aes(color="dSNP per Codon")) +
    #geom_smooth(aes(y=(Number_dSNP/Number_Synonymous)), 
                #method="loess", color="#FC4E07", se=FALSE) +
    #geom_point(aes(y=Number_Tolerated/Num_codons), size=2, shape=23, fill=alpha("#E7B800", 0.8)) +
    #geom_smooth(aes(y=Number_Tolerated/Number_Synonymous), 
                #method="loess", se=FALSE, color="#E7B800") + 
    #geom_point(aes(y=Number_Synonymous/Num_codons), size=2, shape=21, fill=alpha("#00AFBB", 0.8)) +
    theme_minimal() +
    geom_smooth(data=recombination_dataset, aes(x=Mbp, y=cM_Mbp_noOut/scaleNum, color="cM/Mbp"),
                method="loess", se=FALSE, fullrange=TRUE) +
    #scale_y_continuous(name="#Variant/Codon", limit=c(0,NA), 
    #                   sec.axis = sec_axis(~ .*scaleNum, name= "cM/Mbp")) +
    scale_y_continuous(name="", limit=c(0,NA), 
                       sec.axis = sec_axis(~ .*scaleNum, name= "")) +
    coord_cartesian(ylim=c(0,maxY)) +
    labs(x="Mbp", y="", color="Legend") +
    scale_color_manual(values=colors) +
    theme(plot.margin = margin(1,0,0.25,0, "cm")) +
    theme(legend.position = "none") # need this if adding legend in as a ggplot object
  return (plot)
}

Recomb_VarCodonplots <- lapply(names(VariantNums_Chromosome), function(x) {
  Recomb_VarNum_Codon_plot(VariantNums_Chromosome[[x]], 
                       Recomb_ChromCalcs[[x]],
                       0.0015, 7000)})
labels <- paste0("Chromosome ", seq(1,17, by=1))

Recomb_VarCodonplots[[18]] <- as_ggplot(get_legend(dSNP)) # new: add legend as another element inside plot (instead of using "common.legend")

p <- ggarrange(plotlist = Recomb_VarCodonplots, labels=labels, ncol = 4, nrow = 5, label.y = 0.95, 
               hjust = -1, font.label = list(size=12, face="bold"),
               common.legend = FALSE)
               #legend = "right")

annotate_figure(p, left = "Number dSNP/Codon", right = "cM/Mbp", fig.lab = "Deleterious SNPs per Codon")

## 0.002, 7000 to see dSNPs/codon
## 0.02, 700 to see sSNPs & tolerated/codon

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/FigureS2_Recombination_dSNP_codon.pdf")

#width = 8.5, height = 11, units = "in")

#######################################
############### LEGEND ################
#######################################

#### NOTE: can only use the get_legend grob with grid.arrange? (not ggarrange)
### make legend:
ForLegend <- VariantNums_Chromosome$Ha412HOChr01
ForLegend_long <- reshape(ForLegend[,c(2,12,5,9)],
                          direction = "long",
                          varying = c("Number_dSNP", "Number_Tolerated"),
                          v.names = "value",
                          idvar=c("Chromosome", "MiddlePos"),
                          timevar = "Legend", 
                          times=c("Number_dSNP/Num_codons", "cM/Mbp"))

Recomb_Legend <- function(color1, Recombination_color, Y_label) {
  plot <- ggplot(ForLegend_long, aes(x=MiddlePos, y=value)) + 
    geom_smooth(aes(color=Legend), se = FALSE) +
    scale_color_manual(values = c(color1, Recombination_color),
                       labels = c("cM/Mbp", Y_label)) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12))
  return (plot)
}

dSNP <- Recomb_Legend("blue", "#FC4E07", "dSNP/codon")
sSNP <- Recomb_Legend("blue", "#00AFBB", "sSNP/codon")
dSNP_sSNP <- Recomb_Legend("blue", "#D95F02", "dSNP/sSNP")

#######################################
##### RECOMBINATION & sSNP/CODON #######
#######################################

Recomb_sSNP_codon_Plot <- function(variant_dataset, recombination_dataset, maxY, scaleNum) {
  plot <- ggplot(variant_dataset, aes(x=MiddlePos/1000000, y=(Number_Synonymous/Num_codons))) +
    geom_smooth(method="loess", se=FALSE, aes(color="sSNP per Codon")) +
    theme_minimal() +
    geom_smooth(data=recombination_dataset, aes(x=Mbp, y=cM_Mbp_noOut/scaleNum, color="cM/Mbp"),
                method="loess", se=FALSE, fullrange=TRUE) +
    scale_y_continuous(name="", limit=c(0,NA), 
                       sec.axis = sec_axis(~ .*scaleNum, name= "")) +
    coord_cartesian(ylim=c(0,maxY)) +
    labs(x="Mbp", y="", color="Legend") +
    scale_color_manual(values=colors) +
    theme(plot.margin = margin(1,0,0.25,0, "cm")) +
    theme(legend.position = "none") # need this if adding legend in as a ggplot object
  return (plot)
}

Recomb_sSNPCodonplots <- lapply(names(VariantNums_Chromosome), function(x) {
  Recomb_sSNP_codon_Plot(VariantNums_Chromosome[[x]], 
                           Recomb_ChromCalcs[[x]],
                           0.017, 700)})

Recomb_sSNPCodonplots[[18]] <- as_ggplot(get_legend(sSNP))

p2 <- ggarrange(plotlist = Recomb_sSNPCodonplots, labels=labels, ncol = 2, nrow = 9, label.y = 0.95, 
               hjust = -1, font.label = list(size=12, face="bold"),
               common.legend = FALSE)

annotate_figure(p2, left = "Number sSNP/Codon", 
                right = "cM/Mbp", 
                fig.lab = "Synonymous SNPs per Codon"
                )

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/FigureS3_Recombination_sSNP_codon.pdf")

#######################################
##### RECOMBINATION & dSNP/sSNP #######
#######################################

Recomb_dSNP_sSNP_plot <- function(variant_dataset, recombination_dataset, maxY, scaleNum) {
  plot <- ggplot(variant_dataset, aes(x=MiddlePos/1000000, y=(Number_dSNP/Number_Synonymous))) +
    geom_smooth(method="loess", se=FALSE, color="#D95F02") +
    theme_minimal() +
    geom_smooth(data=recombination_dataset, aes(x=Mbp, y=cM_Mbp_noOut/scaleNum, color="cM/Mbp"),
                method="loess", se=FALSE, fullrange=TRUE) +
    scale_y_continuous(name="", limit=c(0,NA), 
                       sec.axis = sec_axis(~ .*scaleNum, name= "")) +
    coord_cartesian(ylim=c(0,maxY)) +
    labs(x="Position (Mbp)", y="", color="Legend") +
    scale_color_manual(values=colors) +
    #theme(plot.margin = margin(1,0,0.25,0, "cm")) +
    theme(legend.position = "none") # need this if adding legend in as a ggplot object
  return (plot)
}

Recomb_dSNPsSNPplots <- lapply(names(VariantNums_Chromosome), function(x) {
  Recomb_dSNP_sSNP_plot(VariantNums_Chromosome[[x]], 
                         Recomb_ChromCalcs[[x]],
                         0.42, 20)})

Recomb_dSNPsSNPplots[[18]] <- as_ggplot(get_legend(dSNP_sSNP))

p3 <- ggarrange(plotlist = Recomb_dSNPsSNPplots, 
                labels=labels, 
                ncol = 2, nrow = 9, 
                #label.y = 0.95, 
                #hjust = -1, 
                #font.label = list(size=12, face="bold"),
                common.legend = FALSE)

annotate_figure(p3, left = "Number dSNP/sSNP", 
                right = "cM/Mbp", 
                #fig.lab = "Number of Deleterious/Synonymous SNPs"
                )

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/FigureS4_dSNP_sSNP_c.pdf",
       width = 4, height = 9, units = "in")


#######################################
####### EXAMPLE W/ CHROMOSOME 7 #######
#######################################

### using R objects from Recombination.R

RecombPlots[[7]]

#ggarrange(Recomb_VarCodonplots[[7]], Recomb_sSNPCodonplots[[7]], Recomb_dSNPsSNPplots[[7]], ncol=1)
#ggarrange(Recomb_VarCodonplots[[1]], Recomb_sSNPCodonplots[[1]], Recomb_dSNPsSNPplots[[1]], ncol=3)

#shapes <- c("dSNP per Codon" = 24, "sSNP per Codon" = 21)

IndChrom <- ggplot(data=VariantNums_Chromosome[[7]], aes(x=MiddlePos/1000000, y=Number_dSNP/Num_codons)) +
  geom_point(size=3, shape=21, aes(fill="dSNP per Codon")) +
  #geom_smooth(method = "loess") +
  geom_point(size=3, shape=21, aes(y=(Number_Synonymous/Num_codons)/10, fill="sSNP per Codon")) +
  scale_y_continuous(name="dSNP/Codon", limit=c(0,NA), 
                     sec.axis = sec_axis(~ .*10, name= "sSNPs/Codon")) +
  labs(x="", color="Legend", shape = "Legend") +
  theme_minimal() +
  theme(legend.title = element_blank(),
    legend.position = "top",
        legend.direction = "horizontal") 

dSNPsSNP_7 <- ggplot(data=VariantNums_Chromosome[[7]], aes(x=MiddlePos/1000000, y=(Number_dSNP/Number_Synonymous))) +
  geom_point(size=3, shape=24, fill="darkgreen") +
  geom_smooth(method="loess", se=FALSE, color="darkgreen") +
  theme_minimal() +
  labs(x="", y="dSNP/sSNP")


### didn't work to change shape in legend manually
#guides(shape=guide_legend(override.aes = list(shape = c(24,21))))
#labs(scale_shape_manual(values = c(24,21)))


ggarrange(IndChrom, dSNPsSNP_7, RecombPlots[[7]], ncol=1, 
          labels=c("A. Variants per codon", "B. Deleterious/Synonymous SNPs", "C. Recombination Rate"), 
          font.label = list(size=12),
          align = "hv")

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Figure3_Chromosome7.pdf")

#######################################
####### BOTH VARIANTS PER CODON #######
#######################################

Var_Per_codon_Plot <- function(variant_dataset, maxY, scaleNum) {
  plot <- ggplot(variant_dataset, aes(x=MiddlePos/1000000, y=(Number_dSNP/Num_codons))) +
    geom_point(size=3, shape=21, aes(fill="dSNP per Codon")) +
    geom_smooth(method="loess", se=FALSE, aes(color="dSNP per Codon")) +
    theme_minimal() +
    geom_point(size=3, shape=21, aes(y=(Number_Synonymous/Num_codons)/scaleNum, fill="sSNP per Codon")) +
    geom_smooth(method="loess", se=FALSE, aes(y=(Number_Synonymous/Num_codons)/scaleNum, color="sSNP per Codon")) +
    scale_y_continuous(name="dSNP/Codon", limit=c(0,maxY), 
                       sec.axis = sec_axis(~ .*scaleNum, name= "sSNPs/Codon")) +
    labs(x="Mbp", y="Number of Variants/Codon", color="Legend") +
    theme(plot.margin = margin(1,0,0.25,0, "cm")) +
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.direction = "horizontal") 
    #theme(legend.position = "none") # need this if adding legend in as a ggplot object
  return (plot)
}


Var_Per_codon_Plot(VariantNums_Chromosome[[7]], 0.00125, 10)

VarperCodonPlots <- lapply(VariantNums_Chromosome, function(x) {
  Var_Per_codon_Plot(x, 0.00125, 10)})

p3 <- ggarrange(plotlist = VarperCodonPlots , labels=labels, ncol = 4, nrow = 5, label.y = 0.95, 
                hjust = -1, font.label = list(size=12, face="bold"),
                common.legend = TRUE)

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/FigureSx_VariantsPerCodon.pdf")

###########SCRATCH BELOW

color="sSNP per Codon"
VariantNums_Chromosome[[7]]



# dSNP/sSNP range 0.07-0.14
Recomb_VarNum_Codon_plot(VariantNums_Chromosome$Ha412HOChr01,
                         Recomb_ChromCalcs$Ha412HOChr01,
                         1, 70)

ggplot(VariantNums_Chromosome$Ha412HOChr01, aes(x=MiddlePos/1000000, y=(Number_dSNP/Num_codons))) +
  geom_point(size=2, shape=24, fill=alpha("#FC4E07", 0.9)) +
  geom_point(aes(y=Number_Tolerated/Num_codons), size=2, shape=23, fill=alpha("#E7B800", 0.8)) +
  geom_point(aes(y=Number_Synonymous/Num_codons), size=2, shape=21, fill=alpha("#00AFBB", 0.8)) +
  theme_minimal() +
  geom_smooth(data=Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp_noOut/700),
              method="loess", se=FALSE, fullrange=TRUE) +
  scale_y_continuous(name="#Variant/Codon", limit=c(0,NA), 
                     sec.axis = sec_axis(~ .*700, name= "cM/Mbp")) +
  coord_cartesian(ylim=c(0,0.015)) +
  xlab("Mbp")


ggplot(VariantNums_Chromosome$Ha412HOChr01, aes(x=MiddlePos/1000000, y=Number_dSNP/Num_codons)) +
  geom_point(size=2, shape=24, fill=alpha("#FC4E07", 0.9)) +
  geom_smooth(method="loess", se=FALSE, color="#FC4E07") +
  geom_smooth(data=Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp_noOut/7000),
              method="loess", se=FALSE, fullrange=TRUE) +
  scale_y_continuous(name="#Variant/Codon", limit=c(0,NA), 
                     sec.axis = sec_axis(~ .*7000, name= "cM/Mbp")) +
  coord_cartesian(ylim=c(0,0.0015))


ggplot(VariantNums_Chromosome$Ha412HOChr01, aes(x=MiddlePos/1000000, y=(Number_dSNP/Number_Synonymous)/Num_codons)) +
  geom_point(size=2, shape=24, fill=alpha("#FC4E07", 0.9)) +
  geom_smooth(data=Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp_noOut/7000),
              method="loess", se=FALSE, fullrange=TRUE) +
  scale_y_continuous(name="#Variant/Codon", limit=c(0,NA), 
                     sec.axis = sec_axis(~ .*7000, name= "cM/Mbp")) +
  coord_cartesian(ylim=c(0,0.0015))

######

colors <- c("dSNP per Codon" = "#FC4E07", "cM/Mbp" = "blue")
Recomb_VarNum_Codon_plot <- function(variant_dataset, recombination_dataset, maxY, scaleNum) {
  plot <- ggplot(variant_dataset, aes(x=MiddlePos/1000000, y=(Number_dSNP/Num_codons))) +
    #geom_point(size=2, shape=24, fill=alpha("#FC4E07", 0.9)) +
    geom_smooth(method="loess", se=FALSE, aes(color="dSNP per Codon")) +
    theme_minimal() +
    geom_smooth(data=recombination_dataset, aes(x=Mbp, y=cM_Mbp_noOut/scaleNum, color="cM/Mbp"),
                method="loess", se=FALSE, fullrange=TRUE) +
    scale_y_continuous(name="", limit=c(0,NA), 
                       sec.axis = sec_axis(~ .*scaleNum, name= "")) +
    coord_cartesian(ylim=c(0,maxY)) +
    labs(x="Mbp", y="", color="Legend") +
    scale_color_manual(values=colors)
  return (plot)
}

Recomb_VarNum_Codon_plot(VariantNums_Chromosome[[1]], 
                           Recomb_ChromCalcs[[1]],
                           0.0015, 7000)

Recomb_VarCodonplots[[1]]


leg_plot <- ggplot(ForLegend_long, aes(x=MiddlePos, y=value))

legP <- leg_plot + geom_smooth(aes(color=Legend), se = FALSE) +
  scale_color_manual(values = c("blue", "#FC4E07")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12))

leg <- get_legend(legP)
ggarrange(Recomb_VarCodonplots[[1]], as_ggplot(leg))
gridExtra::grid.arrange(get_legend(legP), Recomb_VarCodonplots[[1]])

gridExtra::grid.arrange(Recomb_VarCodonplots[[1]], leg)