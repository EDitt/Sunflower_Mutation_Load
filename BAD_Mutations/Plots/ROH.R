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
apply(SNP_ROH, 2, function(x) {
  which(is.na(x))
}) # no na
apply(SNP_ROH, 2, function(x) {
  which(is.infinite(x))
})

SNP_ROH_Numbers <- SNP_ROH[,c(1:5,7,9,11,14)]

# wide to long
SNP_ROH_long <- reshape(SNP_ROH_Numbers,
                        varying=list(names(SNP_ROH_Numbers)[c(5:9)]),
                        v.names = "Value",
                        idvar=c("sample", "Consequence", "NumDerivedHom"),
                        timevar = "Category",
                        times=names(SNP_ROH_Numbers)[c(5:9)],
                        drop=names(SNP_ROH_Numbers)[c(4)],
                        direction = "long"
)


SNP_ROH_long$Proportion <- SNP_ROH_long$Value / SNP_ROH_long$NumDerivedHom

SNP_ROH_long$Category <- factor(SNP_ROH_long$Category,
                                levels=c("Not_in_ROH", "Num_Small", "Num_Medium", "Num_Large", "Tot_HomDer_inROH"))
levels(SNP_ROH_long$Consequence)

SNP_ROH_long$Consequence <- factor(SNP_ROH_long$Consequence,
                                   levels=c("AllDel", "Tolerated", "SynonymousNodups"))

######################################
# BOXPLOT BY CONSEQUENCE AND/OR SIZE #
######################################

# for small, medium, and large roh
ggplot(data=SNP_ROH_long[which(SNP_ROH_long$Category!="Tot_HomDer_inROH"),],
       #aes(x=Category, y=log(Proportion))) +
      aes(x=Category, y=Proportion)) +
  geom_boxplot(aes(color=Consequence), notch = TRUE) +
  ylim(0,0.11)

# across all roh size categories
#SNP_ROH_all <- subset(SNP_ROH_long, Category=="Tot_HomDer_inROH" | Category=="Not_in_ROH")
SNP_ROH_only <- subset(SNP_ROH_long, Category=="Tot_HomDer_inROH")
ggplot(data=SNP_ROH_only,
       aes(x=Consequence, y=Proportion)) +
  geom_boxplot(aes(fill=Consequence), notch = TRUE) +
  #geom_violin(aes(fill=Consequence)) +
  scale_fill_manual(values= c(Col_Deleterious, 
                              Col_Tolerated,
                              Col_Synonymous),
                    labels=c("Deleterious", "Tolerated", "Synonymous"),
                    name = "Variant Class") +
  theme_minimal() +
  xlab("") + ylab("Proportion of SNPs in ROH") +
  scale_x_discrete(labels=c("Deleterious","Tolerated", "Synonymous")) +
  theme(legend.position = "none")
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/SNPProp_ROH_boxplot.pdf",
       width = 8, height = 10)

######################################
########### dSNP/sSNP RATIO ##########
######################################


### need to redo the graph below -- what is SNP_ROH_all???

# look at ratio of dSNP/sSNP 
SNP_ROH_wide <- reshape(SNP_ROH_Numbers[,c(1:2,8,9)],
                        idvar=c("sample"),
                        timevar=c("Consequence"),
                        direction="wide")
SNP_ROH_wide$dSNP_sSNP_inROH <- SNP_ROH_wide$Tot_HomDer_inROH.AllDel / SNP_ROH_wide$Tot_HomDer_inROH.SynonymousNodups
SNP_ROH_wide$dSNP_sSNP_outROH <- SNP_ROH_wide$Not_in_ROH.AllDel / SNP_ROH_wide$Not_in_ROH.SynonymousNodups

SNP_ROH_ratio <- reshape(SNP_ROH_wide[,c(1,8,9)],
                        varying=list(names(SNP_ROH_wide)[c(8:9)]),
                        v.names="Value",
                        idvar=c("sample"),
                        timevar=c("Category"),
                        times=names(SNP_ROH_wide)[c(8:9)],
                        direction="long")
SNP_ROH_ratio$Value[is.na(SNP_ROH_ratio$Value)] <- 0

ggplot(data=SNP_ROH_ratio, aes(x=Category, y=Value)) +
  #geom_violin() +
  geom_boxplot(notch=TRUE) +
  ylim(0,0.15)

# violin plot is the best data for this
ggplot(data=SNP_ROH_ratio, aes(x=Category, y=log(Value))) +
  geom_violin(fill="grey") +
  theme_minimal() +
  ylab("dSNP/sSNP (log-scaled)") + xlab("") +
  scale_x_discrete(labels=c("In ROH","Outside of ROH"))
# ratio of dSNP/sSNP in ROH=1 for PPN018

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/dSNP_sSNP_ROH_violin.pdf",
       width = 8, height = 10)


# differences across germplasm?
SNP_ROH_wgroups <- merge(SNP_ROH_ratio, key,
                         by.x="sample", by.y = "SequenceName")

ggplot(data=SNP_ROH_wgroups, aes(x=group, y=Value)) +
  geom_boxplot(aes(color=Category))

####################################################
### BOXPLOT W/ AVE dSNPs/codon IN AND OUT OF ROH ###
####################################################

SNP_ROH_NumPerCodon <- SNP_ROH[,c(1,2,13,16)]

SNP_ROH_NumPerCodon_long <- reshape(SNP_ROH_NumPerCodon,
                                    varying=list(names(SNP_ROH_NumPerCodon)[c(3:4)]),
                                    v.names="Value",
                                    idvar=c("sample", "Consequence"),
                                    timevar = "Category",
                                    times=names(SNP_ROH_NumPerCodon)[c(3:4)],
                                    direction = "long"
)

# add constant to 0 values?
min(SNP_ROH_NumPerCodon_long[which(SNP_ROH_NumPerCodon_long$Value >0), "Value"]) # 0.000015
length(SNP_ROH_NumPerCodon_long[which(SNP_ROH_NumPerCodon_long$Value == 0), "Value"]) # 32
##
ggplot(data=SNP_ROH_NumPerCodon_long, aes(x=Category, y=log(Value + 0.000001))) +
  geom_boxplot(aes(color=Consequence)) # don't think I can really add a constant here

levels(SNP_ROH_NumPerCodon_long$Consequence)
SNP_ROH_NumPerCodon_long$Consequence <- factor(SNP_ROH_NumPerCodon_long$Consequence,
                                               levels=c("AllDel",
                                                        "Tolerated",
                                                        "SynonymousNodups"))
# this is one possibility for main manuscript
ggplot(data=SNP_ROH_NumPerCodon_long, aes(x=Category, y=log(Value))) +
  geom_violin(aes(fill=Consequence)) +
  #geom_boxplot(aes(fill=Consequence), notch=TRUE) +
  theme_minimal() +
  xlab("") + ylab("Number of Homozygous Derived per Codon (log-scaled)") +
  scale_fill_manual(values= c(Col_Deleterious, 
                              Col_Tolerated,
                              Col_Synonymous),
                    labels=c("Deleterious", "Tolerated", "Synonymous"),
                    name = "Variant Class") +
  scale_x_discrete(labels=c("In ROH","Outside ROH"))
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/SNPNum_Consequence_ROH_violin.pdf",
       width = 8, height = 10)

ggplot(data=SNP_ROH_NumPerCodon_long, aes(x=Category, y=Value)) +
  geom_boxplot() +
  facet_wrap(vars(Consequence), scales="free_y")

ggplot(data=SNP_ROH_NumPerCodon_long, aes(x=Consequence, y=log(Value))) +
  geom_boxplot() +
  facet_wrap(vars(Category), scales="free_y")

# differences across germplasm?
SNP_NumPerCodon_wgroups <- merge(SNP_ROH_NumPerCodon_long, key,
                         by.x="sample", by.y = "SequenceName")

ggplot(data=SNP_NumPerCodon_wgroups[which(SNP_NumPerCodon_wgroups$group!="introgressed"),], 
       aes(x=group, y=log(Value))) +
  geom_violin(aes(fill=Category)) +
  #geom_boxplot(aes(fill=Category)) +
  facet_wrap(vars(Consequence), scales="free_y")

#############################################################
### SCATTERPLOT WITH SUM OF ROH & dSNP/sSNP for ROH SNPs ###
############################################################

SNP_ROH_size <- merge(SNP_ROH, IndROH_INFO,
                      by.x = "sample",
                      by.y = "SequenceName",
                      all.x = TRUE)

SNP_ROH_size$Proportion_inROH <- SNP_ROH_size$Tot_HomDer_inROH / SNP_ROH_size$NumDerivedHom 
SNP_ROH_size$Mb_roh <- SNP_ROH_size$KB/1000

ggplot(data=SNP_ROH_size[which(SNP_ROH_size$Consequence!="Tolerated" &
                                 SNP_ROH_size$group!="introgressed"),], aes(x=KB, y=Proportion_inROH)) +
  #geom_point(aes(color=Groups, shape=Consequence), alpha=0.9, size=3) +
  geom_point(aes(fill=group, shape=Consequence), alpha=0.9, size=3) +
  scale_shape_manual("Consequence", values = c(21,24)) +
  scale_fill_manual("Germplasm Group", values=c(landrace_OPV, landrace_OPV,
                             HA_NonOil, HA_Oil)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
                    #, name="Germplasm Group", labels=c()) +
  geom_smooth(method="lm", aes(linetype=Consequence, color=Consequence), se=FALSE) +
  scale_color_manual(values= c(Col_Deleterious, Col_Synonymous),
                     name = "Variant Class", guide = "none") +
  theme_minimal()

# or just make shape different among germplasm groups?
#ggplot(data=SNP_ROH_size[which(SNP_ROH_size$Consequence!="Tolerated" &
#                                 SNP_ROH_size$group!="introgressed"),], aes(x=KB, y=Proportion_inROH)) +
ggplot(data=SNP_ROH_size[which(SNP_ROH_size$group!="introgressed"),], aes(x=Mb_roh, y=Proportion_inROH)) +
  geom_point(aes(fill=Consequence, shape=group), alpha=0.75, size=3) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  scale_shape_manual("Germplasm Group", values = c(24,23,22,21)) +
  scale_fill_manual("Variant Class", values=c(Col_Deleterious, Col_Synonymous, Col_Tolerated),
                    labels=c("Deleterious", "Synonymous", "Tolerated")) +
  geom_smooth(method="lm", aes(linetype=Consequence, color=Consequence), se=FALSE, show.legend = FALSE) +
  scale_color_manual(name = "Variant Class", 
                     values= c(Col_Deleterious, Col_Synonymous, Col_Tolerated),
                     labels=c("Deleterious", "Synonymous", "Tolerated")) +
  theme_minimal() +
  ylab("Fraction of derived homozygotes in a ROH") +
  xlab("Sum ROH Length (Mb)")
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/ROH_scatterplot.pdf",
       width = 10, height = 10)

ggplot(data=SNP_ROH_size[which(SNP_ROH_size$group!="introgressed"),], aes(x=Mb_roh, y=Proportion_inROH)) +
  geom_smooth(method="lm", aes(linetype=Consequence, color=Consequence), se=FALSE)

######################################################
### SCATTERPLOT WITH dSNP for ROH and non-ROH SNPs ###
######################################################

ggplot(data = SNP_ROH_size[which(SNP_ROH_size$Consequence=="AllDel"),], 
       aes(x=Mb_roh, y=Tot_HomDer_inROH)) +
  geom_point(color="red") + geom_smooth(method="lm", color="red") +
  geom_point(aes(x=Mb_roh, y=Not_in_ROH), color="blue") +
  geom_smooth(aes(x=Mb_roh, y=Not_in_ROH), method="lm", color="blue")

# total dSNPs and roh?
ggplot(data = SNP_ROH_size[which(SNP_ROH_size$Consequence=="AllDel"),], 
       aes(x=Mb_roh, y=NumDerivedHom)) +
  geom_point(color="red") +  geom_smooth(se=TRUE)

# for all consequences?
ggplot(data = SNP_ROH_size,
       aes(x=Mb_roh, y=NumDerivedHom)) +
  geom_point(aes(color=Consequence)) +  
  facet_wrap(vars(Consequence), scales="free_y") +
  geom_smooth(se=TRUE)

# total dSNPs/codon
ggplot(data = SNP_ROH_size[which(SNP_ROH_size$Consequence=="AllDel"),], 
       aes(x=Mb_roh, y=Tot_HomDer_perCodon)) +
  geom_point(color="red") +  geom_smooth(se=TRUE)

ggplot(data = SNP_ROH_size, 
       aes(x=Mb_roh, y=Tot_HomDer_perCodon)) +
  geom_point(aes(color=Consequence)) +  geom_smooth(se=TRUE) +
  facet_wrap(vars(Consequence), scales="free_y")

#########################################################
### BOXPLOT SHOWING TOTAL LENGTH OF ROH OF DIFF SIZES ###
#########################################################




######################################
# N.) Differentiation in proportion of ROH between HA & RHA on chromosome 10

Chr10_ROH <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/Chr10/het_group_summary.txt",
                                  header=TRUE)

Chr10_ROH$HA_prop <- Chr10_ROH$HA_num / 145
Chr10_ROH$RHA_prop <- Chr10_ROH$RHA_num / 112

plot(Chr10_ROH$HA_prop ~ Chr10_ROH$RHA_prop)
