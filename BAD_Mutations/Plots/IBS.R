### scatterplot of IBS values for different crosses (within versus between heterotic groups)

source("BAD_Mutations/Variant_analyses/Functions.R")

library(ggplot2)
library(car)
library(lme4)

######################################
######### SET UP DATAFRAME ###########
######################################

IBS <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/heterotic_groups_new/IBS_syndel.txt",
                  header=T, sep = "\t")

IBS[which(as.character(IBS$Genotype2)==as.character(IBS$Genotype1)),] # check

levels(IBS$Cross)
IBS$Cross <- factor(IBS$Cross, levels = c("Within_Group", "Between_Group"))

######################################
############## MODEL #################
######################################

Mod <- lm(IBS.Deleterious ~ IBS.Synonymous + Cross +
            IBS.Synonymous:Cross,
           data = IBS)
Anova(Mod)
hist(resid(Mod))
summary(Mod)
drop1(Mod, test = "F") # interaction significant, p<0.0001
coef(summary(Mod))

# with random effects?
Mod2 <- lmer(IBS.Deleterious ~ IBS.Synonymous + Cross +
            IBS.Synonymous:Cross +
              (1|Genotype1),
          data = IBS)
Anova(Mod2)
hist(resid(Mod2))
summary(Mod2)
drop1(Mod2, test = "Chisq") # interaction significant, p<0.0001
coef(summary(Mod2))

######################################
######### OVERALL PATTERN ############
######################################

ggplot(data=IBS, aes(x=IBS.Synonymous, y=IBS.Deleterious)) +
  geom_point(aes(color=Cross), alpha=0.2) +
  geom_smooth(method=lm, aes(color=Cross)) +
  geom_abline(slope=1, intercept=0, linetype="dotted") +
  xlim(0.75,1) + ylim(0.75, 1)

######################################
#### DIFFERENCES IN IBS x CROSS ######
######################################

boxplot(IBS$IBS.Deleterious ~ IBS$Cross)
boxplot(IBS$IBS.Synonymous ~ IBS$Cross)

# is there a difference after controlling for IBS at synonymous loci?

Mod_basic <- lm(IBS.Deleterious ~ IBS.Synonymous,
          data = IBS)
boxplot(resid(Mod_basic) ~ IBS$Cross)

######################################
#### USE AVERAGE LINE DISTANCE #######
######################################

MeanDist <- function(df, genotype){
  df_sub <- subset(df, Genotype1==genotype | Genotype2==genotype)
  MeanDel <- aggregate(df_sub$IBS.Deleterious, by=list(df_sub$Cross), mean)
  colnames(MeanDel)[2] <- "MeanIBS_Deleterious"
  MeanSyn <- aggregate(df_sub$IBS.Synonymous, by=list(df_sub$Cross), mean)
  colnames(MeanSyn)[2] <- "MeanIBS_Synonymous"
  All <- merge(MeanDel, MeanSyn, by="Group.1")
  All$Genotype <- genotype
  return(All)
}


geno1 <- unique(as.character(IBS$Genotype1))
All_Geno <- c(geno1, 
              as.character(unique(IBS[which(!IBS$Genotype2 %in% geno1), "Genotype2"])))

MeanGenoDist <- lapply(All_Geno, function(x) {MeanDist(IBS, x)})
MeanGenoDF <- do.call(rbind, MeanGenoDist)
colnames(MeanGenoDF)[1] <- "Cross"

######################################
#### MODEL AVERAGE DIST PER LINE #####
######################################

Mod3 <- lmer(MeanIBS_Deleterious ~ MeanIBS_Synonymous + Cross +
               MeanIBS_Synonymous:Cross +
               (1|Genotype),
             data = MeanGenoDF)
hist(resid(Mod3))

summary(Mod3)
drop1(Mod3, test = "Chi") # n.s.
Anova(Mod3)
boxplot(resid(Mod3) ~ MeanGenoDF$Cross) #slightly lower median for between-group

######################################
##### PLOT AVERAGE DIST PER LINE #####
######################################

ggplot(data=MeanGenoDF, aes(x=MeanIBS_Synonymous, y=MeanIBS_Deleterious)) +
  geom_point(aes(color=Cross), alpha=0.4) +
  geom_smooth(method=lm, aes(color=Cross)) +
  #geom_smooth(method=lm, linetype="dashed", color="black") +
  geom_abline(slope=1, intercept=0, linetype="dotted") +
  theme_minimal() +
  xlim(0.8,0.9) + ylim(0.8, 0.9)


boxplot(MeanGenoDF$MeanIBS_Deleterious ~ MeanGenoDF$Cross, ylab="dSNP distance", ylim=c(0.8, 0.95))
boxplot(MeanGenoDF$MeanIBS_Synonymous ~ MeanGenoDF$Cross, ylab="sSNP distance", ylim=c(0.8, 0.95))

####### SCRATCH BELOW

######################################
############ MAKE PLOT ###############
######################################

ggplot(data=IBS, aes(x=IBS.Synonymous, y=IBS.Deleterious)) +
  geom_point(aes(color=Cross), alpha=0.2) +
  geom_smooth(method=lm, aes(color=Cross)) +
  geom_abline(slope=1, intercept=0, linetype="dotted") +
  geom_abline(slope=coef(summary(Mod2))[,"Estimate"][2], intercept = coef(summary(Mod2))[,"Estimate"][1]) +
  #geom_abline(slope=withingroup, intercept = coef(summary(Mod3))[,"Estimate"][1]) +
  theme_minimal() +
  xlim(0.7,1) + ylim(0.7, 1)



boxplot(IBS$IBS.Deleterious ~ IBS$Cross) ### this looks cool
boxplot(IBS$IBS.Synonymous ~ IBS$Cross) ### similar relationship

Mod3 <- lmer(IBS.Deleterious ~ IBS.Synonymous + Cross +
               IBS.Synonymous:Cross + (1|Genotype1),
             data = IBS)
summary(Mod3)
drop1(Mod3, test = "Chisq") # synonymous x cross p<0.0001
Anova(Mod3) # synonymous x cross more significant than main effect of cross

plot(IBS.Synonymous, IBS.Deleterious, data = IBS)

# within-group slope?
withingroup <- coef(summary(Mod3))[,"Estimate"][2] + coef(summary(Mod3))[,"Estimate"][3] +
  coef(summary(Mod3))[,"Estimate"][4] # 0.7789776

# b/w group intercept = 0.22349
# synonymomus slope = 0.78
# cross within group = 0.009628
# 
coef(summary(Mod3))[,"Estimate"][2]
# intercept and slope are mean across all groups

# plot points separately?
ggplot(data=IBS, aes(x=IBS.Synonymous, y=IBS.Deleterious)) +
  geom_point(data=IBS[which(IBS$Cross=="Between_Group"),], aes(x=IBS.Synonymous, y=IBS.Deleterious), 
             alpha=0.1, color="red") +
  #geom_point(data=IBS[which(IBS$Cross=="Within_Group"),], aes(x=IBS.Synonymous, y=IBS.Deleterious), 
  #           alpha=0.1, color="blue") +
  geom_smooth(method=lm, color="black") +
  theme_minimal()


## duplicate values?
length(IBS$IBS.Deleterious) # 27495
length(unique(IBS$IBS.Deleterious)) # 21720

length(paste0(IBS$Genotype2,"_", IBS$Genotype1)) # 27495
length(unique(paste0(IBS$Genotype2,"_", IBS$Genotype1))) #27495

# average distance for each line to all other HA/RHA?
test_geno <- "PPN001"

IBS_sub <- subset(IBS, Genotype1==test_geno | Genotype2==test_geno)
meandDel1 <- aggregate(IBS_sub$IBS.Deleterious, by=list(IBS_sub$Cross), mean)
colnames(meandDel1)[2] <- "MeanIBS_Deleterious"
MeanSyn1 <- aggregate(IBS_sub$IBS.Synonymous, by=list(IBS_sub$Cross), mean)
colnames(MeanSyn1)[2] <- "MeanIBS_Synonymous"
All <- merge(meandDel1, MeanSyn1, by="Group.1")
All$Genotype <- test_geno




test <- MeanDist(IBS, "PPN001")

genotest <- unique(c(unique(as.character(IBS$Genotype1))), unique(as.character(IBS$Genotype2)))


which(!IBS$Genotype2 %in% genotest)

unique(IBS[which(!IBS$Genotype1 %in% IBS$Genotype2), "Genotype1"]) # this works




