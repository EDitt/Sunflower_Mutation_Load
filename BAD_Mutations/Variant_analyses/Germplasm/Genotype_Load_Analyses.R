# Analyses to look at patterns of dSNPs among germplasm
##  - the Distribution among different germplasm groups
##  - relationship with performance

#### continues from "GenotypeLoadNums.md"

####################
###### SETUP #######
####################

library(car)
library(emmeans)
library(ggplot2)
library(reshape2)

setwd("/Users/emilydittmar/Google Drive/Active Projects/DelMutation/Results/Genotype_patterns")

######################################
######### SET UP DATAFRAME ###########
######################################

# see GenotypeLoadNums.md for how this file was generated
dSNP_genotypes <- read.table("Genotype_dSNP_counts.txt", sep = "\t", header=TRUE,
                            stringsAsFactors = FALSE)

# Genotype Information
key <- read.csv("LineKeywINFO.csv", header=T)

dSNP_GenoInfo <- merge(key, dSNP_genotypes, by.x="SequenceName", by.y="sample")
# check merge
length(dSNP_GenoInfo$SequenceName) # 288

# samples have different numbers of called genotypes
### if I want to relativize by number of called genotypes

dSNP_GenoInfo$NumGenotypes <- dSNP_GenoInfo$Total - dSNP_GenoInfo$TotalMissing

#### Note: all numbers correspond to the numbers of *genotypes* 
### except TotNum_dSNPs which refers to number of *alleles* (ie 2*nHom)

str(dSNP_GenoInfo)
levels(dSNP_GenoInfo$heterotic_group)
levels(dSNP_GenoInfo$Oil_NonOil)
levels(dSNP_GenoInfo$type)

# look at nHom_dSNPs relativized by called genotypes
dSNP_GenoInfo$nHom_dSNPs_REL <- dSNP_GenoInfo$nHom_dSNPs / dSNP_GenoInfo$NumGenotypes
hist(dSNP_GenoInfo$nHom_dSNPs_REL)
dSNP_GenoInfo[which(dSNP_GenoInfo$nHom_dSNPs_REL < 0.1),] # Hopi

####################
###### ANOVA #######
####################

# investigate whether categories effect TotNum_dSNPs or nHom_dSNPs
hist(dSNP_GenoInfo$TotNum_dSNPs)

mod1 <- glm(TotNum_dSNPs ~ offset(log(2*NumGenotypes)) +
              heterotic_group + Oil_NonOil+
             heterotic_group:Oil_NonOil,
            family = poisson(link = "log"),
           data=dSNP_GenoInfo)

summary(mod1)
Anova(mod1, test= "F") # all categories are significant
hist(resid(mod1)) # normal

Mod_means <- emmeans(mod1, ~ heterotic_group | Oil_NonOil, type = "response")
Mod_means
pairs(Mod_means)

# heterotic group?
pairs(emmeans(mod1, ~ heterotic_group, type = "response"))

emmip(mod1, ~ heterotic_group | Oil_NonOil, type = "response")

Mean.df <- as.data.frame(Mod_means)
Mean.df$group <- paste0(Mean.df$heterotic_group,"-", Mean.df$Oil_NonOil)
levels(as.factor(Mean.df$group))
Mean.df_sub <- subset(Mean.df, !is.na(rate))
Mean.df_sub$group <- factor(Mean.df_sub$group, levels=c(
  "HA-INRA", "RHA-INRA", "HA-NonOil", "HA-Oil", "RHA-NonOil",
  "RHA-Oil", "introgressed-NonOil", "introgressed-Oil", "landrace-NonOil",
  "OPV-NonOil", "other-NonOil", "other-Oil"
))
levels(as.factor(Mean.df_sub$group))

p <- ggplot(data=Mean.df_sub, mapping=aes(x=group, y=rate))
p + geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=rate-SE, ymax=rate+SE,
                    width=.2))

# just look at type
mod2 <- glm(TotNum_dSNPs ~ offset(log(2*NumGenotypes)) +
              type,
            family = poisson(link = "log"),
            data=dSNP_GenoInfo)

summary(mod2)
Anova(mod2, test= "F")
hist(resid(mod2)) # normal

Mod_means2 <- emmeans(mod2, ~ type, type = "response")
Mod_means2
pairs(Mod_means2)
emmip(mod2, ~ type, type = "response")

# what about # of homozygous deleterious genotypes?
mod3 <- glm(nHom_dSNPs ~ offset(log(NumGenotypes)) +
              type,
            family = poisson(link = "log"),
            data=dSNP_GenoInfo)
emmip(mod3, ~ type, type = "response")


####################
### PERFORMANCE ####
####################

performance <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Nut_Pilot/data files/HELNUT16SAMlsmeans_wide.csv", header=T)
length(performance$SAM) #263

AllData <- merge(performance[,c(2,43,44,53,54,59,60,419)], dSNP_GenoInfo, by.y="SamID", by.x="SAM")
length(AllData$SAM) #262

# long format
Data_long <- melt(AllData[,c(1:8,17:19,21:22,26,27)], id.vars = c("SAM", "type", "heterotic_group",
                                                         "Oil_NonOil", "TotNum_dSNPs",
                                                         "nHom_dSNPs", "NumGenotypes", "nHom_dSNPs_REL") )
Data_long$variable <- strsplit(as.character(Data_long$variable), "_")
Data_long$PerfMeasure <- as.factor(sapply(Data_long$variable, "[", 1))
Data_long$Condition <- as.factor(sapply(Data_long$variable, "[", 2))

### model function
Model_dSNP_Performance <- function (dataframe, performance_metric, dSNP_metric) {
  newdf <- subset(dataframe, PerfMeasure==performance_metric)
  newdSNP_metric <- dataframe[which(Data_long$PerfMeasure==performance_metric),dSNP_metric]
  model <- lm(value ~ newdSNP_metric  + Condition +
                heterotic_group +
                newdSNP_metric:Condition,
              data=newdf)
  return(model)
}

# added RGR diff
Mod1 <- Model_dSNP_Performance(Data_long[which(Data_long$Condition!="diff")], "RGR", "TotNum_dSNPs")
summary(Mod1)
drop1(Mod1, test="Chi") #dSNPs x Condition ns
Anova(Mod1) # only condition and heterotic group are sig.

# homozygous genotypes
Mod2 <- Model_dSNP_Performance(Data_long[which(Data_long$Condition!="diff")], "RGR", "nHom_dSNPs_REL")
summary(Mod2)
drop1(Mod2, test="Chi") #nHomdSNPs x Condition p=0.005
Anova(Mod2) # All sig.

# plant weight?
Mod4 <- Model_dSNP_Performance(Data_long[which(Data_long$Condition!="diff")], "ln.Plant.weight", "nHom_dSNPs_REL")
summary(Mod4)
drop1(Mod4, test="Chi") #ratio nHom_dSNPs x condition is ns
Anova(Mod4) # dSNP metric is marginally significant (p=0.05)


#####################
# PLOT RELATIONSHIP #
#####################

Plot_dSNP_Performance <- function (dataframe, performance_metric, dSNP_metric) {
  newdf <- subset(dataframe, PerfMeasure==performance_metric)
  newdSNP_metric <- dataframe[which(Data_long$PerfMeasure==performance_metric),dSNP_metric]
  p <- ggplot(data=newdf, aes(x=newdSNP_metric, y=value, group=Condition)) +
  geom_point(aes(color=Condition)) +
    geom_smooth(method='lm', aes(color=Condition)) +
    theme_minimal() +
    ylab("Relative growth rate")
  return(p)
}

## adding RGR diff screwed things up
Plot_dSNP_Performance(Data_long[which(Data_long$Condition!="diff")], "RGR", "TotNum_dSNPs")
Plot_dSNP_Performance(Data_long[which(Data_long$Condition!="diff")], "RGR", "nHom_dSNPs")
Plot_dSNP_Performance(Data_long[which(Data_long$Condition!="diff")], "RGR", "nHom_dSNPs_REL")
Plot_dSNP_Performance(Data_long[which(Data_long$Condition!="diff")], "ln.Plant.weight", "nHom_dSNPs_REL")






### unused so far
# stress / control
AllData$StressControl <- AllData$RGR_Low / AllData$RGR_High

plot(AllData$RGR_High ~ AllData$TotNum_dSNPs)
plot(AllData$RGR_Low ~ AllData$TotNum_dSNPs)

## older code (already incorporated)
# relativize by called genotypes?
Data_long$RATIO_nHom_dSNPs <- Data_long$nHom_dSNPs / Data_long$NumGenotypes
Mod3 <- Model_dSNP_Performance(Data_long, "RGR", "RATIO_nHom_dSNPs")
summary(Mod3)
drop1(Mod3, test="Chi") #ratio nHom_dSNPs x condition is significant (p=0.0053)
Anova(Mod3) # ratio nHom_dSNPs sig. (p=0.03)


