### Plotting runs of homozygosity

source("BAD_Mutations/Variant_analyses/Functions.R")

# 1.) "hotspots" of ROH, i.e. stretches of homozygous sequence shared by a large proportion of individuals

# see README.md -> Runs of homozygosity for how data were created

######################################
######### SET UP DATAFRAME ###########
######################################

#GenomeROH <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/ROH_freq_1MbpBins.txt",
GenomeROH <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/LROH_2Mbp_freq_1MbpBins.txt",
                        header=F, 
                        col.names=c("Chromosome", "StartPos", "EndPos",
                                    "Mean_Num_Indiv", "Median_Num_Indiv",
                                   "Num_SNPs"))

GenomeROH$Mbp <- round((GenomeROH$StartPos / 1000000) + 1, digits = 0)

######################################
######## PLOT ACROSS GENOME ##########
######################################

hist(GenomeROH$Mean_Num_Indiv)
hist(GenomeROH$Median_Num_Indiv)

plotFst(GenomeROH, "Median_Num_Indiv", c(0.975, 0.995))
GenomeROH[which(GenomeROH$Median_Num_Indiv > 150),]

plotFst(GenomeROH, "Mean_Num_Indiv", c(0.975, 0.995))


######################################
# 2.) length and number of ROH for all individuals (color-coded by group)

######################################
######### SET UP DATAFRAME ###########
######################################

IndROH <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/Plink_default.hom.indiv",
                        header=TRUE)
IndROH$FID <- as.character(IndROH$FID)
IndROH$IID <- as.character(IndROH$IID)

# fix sequence name
IndROH$SequenceName <- ifelse(IndROH$FID==IndROH$IID,
                             as.character(IndROH$FID),
                             as.character(paste0(IndROH$FID, "_", 
                                                 IndROH$IID)))

# set up groupings
key <- read.csv("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Sunflower_Mutation_Load/BAD_Mutations/LineKeywINFO.csv", 
                header=T)
key$group <- ifelse(key$heterotic_group=="landrace" |
                        key$heterotic_group=="OPV",
                      as.character(key$heterotic_group),
                      as.character(key$Groups))

# combine HA/RHA- keep as Oil vs. NonOil
key$group <- ifelse(key$heterotic_group=="landrace",
                    "landrace",
                    ifelse(key$heterotic_group=="HA" |
                             key$heterotic_group=="RHA",
                           as.character(key$Oil_NonOil),
                           as.character(key$Groups)))

key$group[which(key$group=="landrace_OPV")] <- "OPV"

levels(as.factor(key$group))
key$group <- factor(key$group, levels=c("landrace", "OPV", "introgressed",
                                        "NonOil", "Oil"))

IndROH_INFO <- merge(IndROH[,c(7,4,5,6)],
                        key,
                        by = "SequenceName")

IndROH_INFO$Mb <- IndROH_INFO$KB / 1000

######################################
###### PLOT BY GERMPLASM GROUP #######
######################################

ggplot(data=IndROH_INFO, aes(x=IndROH_INFO$NSEG,
                             y=IndROH_INFO$KB)) +
  geom_point(aes(fill=group, shape=group, color=group), size=3, alpha=0.8) +
  scale_shape_manual(values=c(8, 23, 21, 21, 21), name = "Germplasm Group") +
  scale_fill_manual(limits=c("landrace", "OPV", "NonOil", "Oil"),
                       values= c(landrace_OPV, landrace_OPV, 
                                 HA_NonOil, HA_Oil),
                       name = "Germplasm Group") +
  scale_color_manual(values=c(landrace_OPV, "black", "black", "black", "black")) +
  theme_minimal() +
  xlab("nROH") + ylab("Total Length")


# to color-code:
ggplot(data=IndROH_INFO[which(IndROH_INFO$group!="landrace"),], 
        aes(x=group, y=Mb)) +
  #geom_violin(outlier.colour = NULL, aes(colour = group))
  geom_boxplot(notch = FALSE, outlier.colour = NULL, aes(colour = group)) +
  geom_boxplot(outlier.shape = NA, aes(fill = group)) +  # in order to keep outline of box black but have outlier colors the same
  theme_minimal() +
  ylab("Total ROH Length") +
  scale_fill_manual(values= c(landrace_OPV, 
                              HA_NonOil, HA_Oil, landrace_OPV),
                    name = "Variant Class") +
  scale_color_manual(values= c("black", 
                              HA_NonOil, HA_Oil, landrace_OPV),
                    name = "Variant Class") +
  scale_x_discrete(limits=c("landrace", "OPV", "NonOil", "Oil")) +
  geom_point(data=IndROH_INFO[which(IndROH_INFO$group=="landrace"),], 
             position = "identity",
             aes(x=1, y=Mb, color="black", fill=landrace_OPV),
             shape=c(22,23,24), size = 3)

Burden_boxplot(IndROH_INFO, "group", "Mb", "Total Length of ROH (Mb)", 
               c("landrace", "OPV", "NonOil", "Oil"))
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/ROH_Length_Boxplot.pdf",
       width = 5, height = 7)

Burden_boxplot(IndROH_INFO, "Groups", "Mb", "Total Length of ROH (Mb)", 
               c("landrace", "OPV", 
                 "HA-NonOil", "RHA-NonOil",
                 "HA-Oil", "RHA-Oil"))

######################################
# 3.) how many dSNPs are in ROH?

######################################
######### SET UP DATAFRAME ###########
######################################

SNP_ROH <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/ROH_SNP_info.txt",
                     header=TRUE)

# wide to long
SNP_ROH_long <- reshape(SNP_ROH,
                        varying=list(names(SNP_ROH)[c(5:9)]),
                        v.names = "Value",
                        idvar=c("sample", "Consequence", "NumDerivedHom"),
                        timevar = "Category",
                        times=names(SNP_ROH)[c(5:9)],
                        drop=names(SNP_ROH)[c(4)],
                        direction = "long"
)
SNP_ROH_long$Proportion <- SNP_ROH_long$Value / SNP_ROH_long$NumDerivedHom

######################################
### BOXPLOT BY CONSEQUENCE AND SIZE ###
######################################

ggplot(data=SNP_ROH_long[which(SNP_ROH_long$Category!="Tot_HomDer_inROH"),],
       aes(x=Category, y=Proportion)) +
  geom_boxplot(aes(color=Consequence))


######################################
# N.) Differentiation in proportion of ROH between HA & RHA on chromosome 10

Chr10_ROH <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/Chr10/het_group_summary.txt",
                                  header=TRUE)

Chr10_ROH$HA_prop <- Chr10_ROH$HA_num / 145
Chr10_ROH$RHA_prop <- Chr10_ROH$RHA_num / 112

plot(Chr10_ROH$HA_prop ~ Chr10_ROH$RHA_prop)
