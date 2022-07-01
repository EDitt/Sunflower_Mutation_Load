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
  

######################################
############# FUNCTION ###############
######################################
  
plotFst <- function(dataframe, ColName){
  quantiles <- quantile(dataframe[,ColName], c(0.975, 0.995), na.rm = T)
  dataframe$PointCol <- ifelse(dataframe[,ColName]>quantiles[2],
                                "99%",
                                ifelse(dataframe[,ColName]>quantiles[1],
                                       "95%",
                                       ifelse((as.numeric(gsub("Ha412HOChr", "", dataframe$Chromosome)) %%2)==0,
                                              "even",
                                              "odd")))
  dataframe$PointCol <- factor(dataframe$PointCol, levels = c("99%", "95%", "odd", "even"))
  ChromosomeNames <- as.character(unique(as.numeric(gsub("Ha412HOChr", "", dataframe$Chromosome))))
  names(ChromosomeNames) <- unique(dataframe$Chromosome)
  plot <- ggplot(data=dataframe, aes(x=Mbp, y=dataframe[,ColName])) +
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
    xlab("Chromosome Position") +
    ylab("Fst")
  return(plot)
}

  
######################################
######## PLOT ACROSS GENOME ##########
######################################


plotFst(Fst_HA_RHA, "Weighted_Fst")
plotFst(subset(Fst_df, Variant_type=="Oil_NonOil"), "Weighted_Fst")


######################################
########### TROUBLESHOOT #############
######################################
Fst_df$Mbp <- round((Fst_df$StartPos / 1000000) + 1, digits = 0) # need to add 5 otherwise will not graph 0 class
ggplot(data=Fst_df, aes(x=Mbp, y=Mean_Fst)) +
  geom_line(aes(color=Variant_type)) +
  facet_grid(~Chromosome, scales = "free_x", space = "free_x")
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