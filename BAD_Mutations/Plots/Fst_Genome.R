### plotting Fst across the genome
### Fst calculated in 1 Mbp windows across the genome

source("BAD_Mutations/Variant_analyses/Functions.R")

library(ggplot2)
library(tidyverse)

# weighted Fst: derived from the ratio of the separate per-locus sums of numerator and denominator values
    # sum(a) / sum(a+b) , i.e. ratio of averages <- recommended by Bhatia et al. 2013
      # weighted average found to provide low bias and variance (see Weir and Hill 2002)
# unweighted "average" Fst: average per-locus values
    # average(a/(a+b)) , i.e. average of ratios

######################################
######### SET UP DATAFRAME ###########
######################################

Fst_df <- ImportFilesAsDf("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicPatterns/Fst_data",
                             ".windowed.weir.fst", "", 
                          c("Chromosome", "StartPos", "EndPos",
                            "N_Variants", "Weighted_Fst", "Mean_Fst"),
                          TRUE)

Fst_df$Variant_type[grep("HA_RHA", Fst_df$Variant_type)] <- "HA_RHA"
Fst_df$Variant_type[grep("Oil_NonOil", Fst_df$Variant_type)] <- "Oil_NonOil"

Fst_df$Mbp <- round((Fst_df$StartPos / 1000000) + 1, digits = 0)  # need to add 1 otherwise will not graph 0 class

Fst_df[which(Fst_df$Weighted_Fst > quantile(Fst_df$Weighted_Fst, 0.995)),] # 0.42
# on chrom 10: 7, 10-25, 30-33, 38-39, 45-46, 90,109

Fst_df[which(Fst_df$Weighted_Fst > quantile(Fst_df$Weighted_Fst, 0.975)),] # 0.268
# on chrom 10: 1-2,6-40, 42-48, 52, 69, 71-82, 89-112, 117-120, 122-124, 
#               126-128, 130, 136, 140, 145, 149-150, 156-157, 162, 169
# on chrom 13: 159-166, 174
######################################
######## PLOT ACROSS GENOME ##########
######################################

ChromosomeNames <- as.character(unique(as.numeric(gsub("Ha412HOChr", "", Fst_HA_RHA$Chromosome))))
names(ChromosomeNames) <- unique(Fst_HA_RHA$Chromosome)

plotFst(subset(Fst_df, Variant_type=="HA_RHA"), "Weighted_Fst", c(0.975, 0.995))
Fst_df[which(Fst_df$Weighted_Fst > 0.4 & Fst_df$Variant_type=="HA_RHA"),]

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Fst_HA_RHA.pdf",
       width = 10, height = 10)

plotFst(subset(Fst_df, Variant_type=="Oil_NonOil"), "Weighted_Fst", c(0.975, 0.995))
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Fst_Oil_NonOil.pdf",
       width = 10, height = 10)

######################################
######### PLOT BOTH AS GRID ###########
######################################

ggplot(data=Fst_df, aes(x=Mbp, y=Weighted_Fst)) +
  geom_line(aes(color=Variant_type)) +
  facet_wrap(~Chromosome, scales = "free_x")


######################################
######### IMPORT AND GRAPH PI #########
######################################

Pi <- ImportFilesAsList("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicPatterns/Fst_data",
                        ".windowed.pi", "_Mbp1_", TRUE)
Pi_HA <- Pi[c(grep("^HAH", names(Pi)))]
names(Pi_HA) <- gsub("HA", "", names(Pi_HA))

Pi_HA <- lapply(Pi_HA, function(x) {
  x$group <- "HA"; return(x)
})

Pi_RHA <- Pi[c(grep("^RHA", names(Pi)))]
names(Pi_RHA) <- gsub("RHA", "", names(Pi_RHA))

Pi_RHA <- lapply(Pi_RHA, function(x) {
  x$group <- "RHA"; return(x)
}) 

Pi_all <- rbind(do.call("rbind", Pi_HA),
             do.call("rbind", Pi_RHA))

Pi_all$Mbp <- round((Pi_all$BIN_START / 1000000) + 1, digits = 0) 
Pi_all <- Pi_all[order(Pi_all$CHROM, Pi_all$Mbp),]


ggplot(data=Pi_all, aes(x=Mbp, y=PI)) +
  geom_smooth(aes(color=group)) +
       #geom_line(aes(color=group)) +
  facet_wrap(~CHROM, scales = "free_x")


# for wide format
Pi_chrom <- lapply(names(Pi_HA), function(x) {
  merge(Pi_HA[[x]], Pi_RHA[[x]], by=c("CHROM", "BIN_START", "BIN_END"), suffixes=c("_HA", "_RHA"))
})

Pi_all_wide <- do.call("rbind", Pi_chrom)
Pi_all_wide$Mbp <- round((Pi_all_wide$BIN_START / 1000000) + 1, digits = 0) 
Pi_all_wide$PiHA_PiRHA <- Pi_all_wide$PI_HA / Pi_all_wide$PI_RHA
Pi_all_wide$PiRHA_PiHA <- Pi_all_wide$PI_RHA / Pi_all_wide$PI_HA

ggplot(data=Pi_all_wide, aes(x=Mbp, y=PiRHA_PiHA)) +
  #geom_line() +
  geom_smooth() +
  facet_wrap(~CHROM, scales = "free_x")

ggplot(data=Pi_all_wide, aes(x=Mbp, y=PiRHA_PiHA)) +
  #geom_line() +
  geom_smooth() +
  facet_grid(~CHROM, scales = "free_x")

######################################
#### IMPORT AND GRAPH TAJIMA'S D #####
######################################

HetGroups_Dem <- function(DIR, suffix, prefix, headerTF) {
  Data <- ImportFilesAsList(DIR, suffix, prefix, headerTF)
  Data_HA <- Data[c(grep("^HAH", names(Data)))]
  names(Data_HA) <- gsub("HA", "", names(Data_HA))
  Data_HA <- lapply(Data_HA, function(x) {
    x$group <- "HA"; return(x)
  })
  Data_RHA <- Data[c(grep("^RHA", names(Data)))]
  names(Data_RHA) <- gsub("RHA", "", names(Data_RHA))
  Data_RHA <- lapply(Data_RHA, function(x) {
    x$group <- "RHA"; return(x)
  })
  Data_all <- rbind(do.call("rbind", Data_HA),
                  do.call("rbind", Data_RHA))
  Data_all$Mbp <- round((Data_all$BIN_START / 1000000) + 1, digits = 0) 
  Data_all <- Data_all[order(Data_all$CHROM, Data_all$Mbp),]
}


Pi_test <- HetGroups_Dem("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicPatterns/Fst_data",
                        ".windowed.pi", "_Mbp1_", TRUE)

ggplot(data=Pi_test, aes(x=Mbp, y=PI)) +
  geom_smooth(aes(color=group)) +
  #geom_line(aes(color=group)) +
  facet_wrap(~CHROM, scales = "free_x")
 

TajD <- HetGroups_Dem("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicPatterns/Fst_data",
                         ".Tajima.D", "_Mbp1_", TRUE)

ggplot(data=TajD, aes(x=Mbp, y=TajimaD)) +
  geom_smooth(aes(color=group)) +
  #geom_line(aes(color=group)) +
  facet_wrap(~CHROM, scales = "free_x")

ggplot(data=TajD, aes(x=Mbp, y=TajimaD)) +
  geom_smooth(aes(color=group)) +
  #geom_line(aes(color=group)) +
  facet_grid(~CHROM, scales = "free_x")


######################################
########### TROUBLESHOOT #############
######################################




Pi_all_long <- reshape(Pi_all,
                       varying=c("PI_HA"))


Pi_df <- ImportFilesAsDf("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicPatterns/Fst_data",
                         ".windowed.pi", "Pi_1Mbp__", 
                         c("Chromosome", "StartPos", "EndPos",
                           "N_Variants", "Pi"),
                         TRUE)

Pi_df$Mbp <- round((Pi_df$StartPos / 1000000) + 1, digits = 0) 
Pi_quant <- quantile(Pi_df$Pi, probs=c(0.05,0.95), na.rm = T)

ggplot(data=Pi_df, aes(x=Mbp, y=Pi)) +
  geom_line() +
  geom_hline(yintercept = c(Pi_quant)) +
  facet_wrap(~Chromosome, scales = "free_x")

plotFst(Pi_df, "Pi", c(0.05,0.95))

###
ggplot(data=Fst_df, aes(x=Mbp, y=Mean_Fst)) +
  geom_line(aes(color=Variant_type)) +
  #facet_grid(~Chromosome, scales = "free_x", space = "free_x")
facet_wrap(~Chromosome, scales = "free_x")

ggplot(data=Fst_df, aes(x=Mbp, y=Weighted_Fst)) +
  geom_line(aes(color=Variant_type)) +
  facet_grid(~Chromosome, scales = "free_x", space = "free_x") +
  hline(Het_threshold)


# just HA/RHA
ggplot(data=Fst_HA_RHA, aes(x=Mbp, y=Weighted_Fst)) +
  #geom_line() +
  geom_point(aes(color=PointCol)) +
  facet_grid(~Chromosome, 
             scales = "free_x", space = "free_x", switch = "x",
             labeller = labeller(Chromosome=ChromosomeNames)) +
  theme_minimal() +
  scale_color_manual(
    values= c("tomato1", "darkgoldenrod1",
              "lightskyblue3", "grey70"),
    name = "", labels = c("99% quantile", 
                          "95% quantile", "", "")) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.spacing.x = unit(0, "null"),
        axis.text.x = element_blank()) +
  xlab("Chromosome Position")


ChromosomeNames <- as.character(unique(as.numeric(gsub("Ha412HOChr", "", Fst_HA_RHA$Chromosome))))
names(ChromosomeNames) <- unique(Fst_HA_RHA$Chromosome)

#  geom_hline(yintercept = c(0.2, 0.4, 0.6, 0.8), color = "grey", linetype = "solid")


#  scale_shape_manual(values = c(8, 19, 19, 19))

######################################
############# FST STATS ##############
######################################

Fst_HA_RHA <- subset(Fst_df, Variant_type=="HA_RHA")
mean(Fst_HA_RHA$Mean_Fst) # 0.0274483
mean(Fst_HA_RHA$Weighted_Fst) # 0.06
quantiles <- quantile(Fst_HA_RHA$Weighted_Fst, c(0.975, 0.995), na.rm = T)
Het_threshold <- quantile(Fst_HA_RHA$Mean_Fst, 0.975, na.rm = T) # 95% quantile 0.153
#Het_threshold <- quantile(Fst_HA_RHA$Mean_Fst, 0.995, na.rm = T) # 99% quantile # 0.27

Fst_HA_RHA[which(Fst_HA_RHA$Mean_Fst>=Het_threshold),]

Fst_HA_RHA$Quantile <- ifelse(Fst_HA_RHA$Weighted_Fst>quantiles[2],
                              "99%",
                              ifelse(Fst_HA_RHA$Weighted_Fst>quantiles[1],
                                     "95%",
                                     "background"))

levels(as.factor(Fst_HA_RHA$Quantile))

Fst_HA_RHA$ChrNum <- as.numeric(gsub("Ha412HOChr", "", Fst_HA_RHA$Chromosome))

Fst_HA_RHA$ChrCol <- ifelse((as.numeric(gsub("Ha412HOChr", "", Fst_HA_RHA$Chromosome)) %%2)==0,
                            "even",
                            "odd")


aggregate(Fst_HA_RHA$StartPos, by=list(Fst_HA_RHA$Chromosome, Fst_HA_RHA$ChrCol), length)


# assign colors to points
Fst_HA_RHA$PointCol <- ifelse(Fst_HA_RHA$Weighted_Fst>quantiles[2],
                              "99%",
                              ifelse(Fst_HA_RHA$Weighted_Fst>quantiles[1],
                                     "95%",
                                     ifelse((as.numeric(gsub("Ha412HOChr", "", Fst_HA_RHA$Chromosome)) %%2)==0,
                                            "even",
                                            "odd")))

aggregate(Fst_HA_RHA$StartPos, by=list(Fst_HA_RHA$Chromosome, Fst_HA_RHA$PointCol), length)

levels(as.factor(Fst_HA_RHA$PointCol))
Fst_HA_RHA$PointCol <- factor(Fst_HA_RHA$PointCol, levels = c("99%", "95%", "odd", "even"))