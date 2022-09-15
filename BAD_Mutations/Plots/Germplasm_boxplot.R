### 1.) looking at average number of dSNPs per genotype and graphing additive, dominant, recessive load per germplasm groups
## 2.) derived dSNPs/sSNPs for groups of germplasm

source("BAD_Mutations/Variant_analyses/Functions.R")

# formerly Genotype_Load_Analyses.R

library(car)
library(emmeans)

#####################################
######### SET UP GROUPINGS ##########
#####################################

key <- read.csv("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Sunflower_Mutation_Load/BAD_Mutations/LineKeywINFO.csv", 
                header=T)

# separate landrace & OPV

key$Groups2 <- ifelse(key$heterotic_group=="landrace" |
                        key$heterotic_group=="OPV",
                      as.character(key$heterotic_group),
                      as.character(key$Groups))
levels(as.factor(key$Groups2))
key$Groups2 <- factor(key$Groups2, levels=c("landrace", "OPV", "introgressed",
                                            "HA-NonOil", "RHA-NonOil",
                                            "HA-Oil", "RHA-Oil"))

# combine HA/RHA- keep as Oil vs. NonOil
key$group <- ifelse(key$heterotic_group=="landrace",
                    "landrace",
                    ifelse(key$heterotic_group=="HA" |
                             key$heterotic_group=="RHA",
                           as.character(key$Oil_NonOil),
                           as.character(key$Groups)))

key$group[which(key$group=="landrace_OPV")] <- "OPV"


aggregate(key$sort, by=list(key$group, key$Groups), length) # check
aggregate(key$sort, by=list(key$group), length) # sample sizes

levels(as.factor(key$group))
key$group <- factor(key$group, levels=c("landrace", "OPV", "introgressed",
                                        "NonOil", "Oil"))

###################################### 1.) All dSNP data

######################################
######### READ IN DATAFRAME ##########
######################################

dSNP_nums <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Genotype_patterns/All_dSNP_stats.txt",
                        sep = "\t", header=T)

######################################
######## AVERAGE dSNP per IND ########
######################################

dSNP_nums$nDeleterious_Total <- 2*dSNP_nums$nDeleterious_Hom +
  dSNP_nums$nHet

hist(dSNP_nums$nDeleterious_Total)

mean(dSNP_nums$nDeleterious_Total) # 23,329.93
median(dSNP_nums$nDeleterious_Total) # 23,771.5

mean(dSNP_nums$nDeleterious_Hom) # 11,055.08
mean(dSNP_nums$nHet) # 1,219.774

# mean number of sites with a deleterious allele
mean(dSNP_nums$nDeleterious_Hom + dSNP_nums$nHet) # 12274.85
mean((dSNP_nums$nDeleterious_Hom + dSNP_nums$nHet)/dSNP_nums$CalledGenotypes) # 0.0003694593, 0.0369 %

# proportion of called genotypes
mean (dSNP_nums$nDeleterious_Total / (2*dSNP_nums$CalledGenotypes)) # 0.00035, 0.035 %

mean(dSNP_nums$nHet / dSNP_nums$CalledGenotypes) #  0.0000367, 0.00367%
mean(dSNP_nums$nDeleterious_Hom / dSNP_nums$CalledGenotypes) # 0.0003327 0.03327%

mean(dSNP_nums$nHet / dSNP_nums$nDeleterious_Hom) # 0.119

#####################################
######### MUTATIONAL BURDEN #########
#####################################

### count of derived deleterious alleles carried by an individual divided by the total number of scored (nonmissing) alleles - Valluru et al. 2019
### additive burden: sum of derived alleles at deleterious sites, (see Lozano et al. )
#####   under dominant model, each site with at least 1 del allele counted, 
#####   under additive model, the number of alleles that were deleterious
#####   under recessive model, only sites that were homozygous for deleterious variants were counted
#####   also used normalized additive model which relativized by the number of available deleterious loci (to account for differences in call rates)

# will relative all metrics by the total number of scored alleles or genotypes:

##### additive: count of total number of deleterious alleles relatived by # of scored alleles
dSNP_nums$Rel_additive <- dSNP_nums$nDeleterious_Total / (2*dSNP_nums$CalledGenotypes)
hist(dSNP_nums$Rel_additive)

##### dominant: count of sites containing at least 1 deleterious allele relatived by # of called genotypes
dSNP_nums$Rel_dominant <- (dSNP_nums$nDeleterious_Hom + dSNP_nums$nHet) / dSNP_nums$CalledGenotypes

##### recessive: count of homozygous deleterious sites relatived by # of called genotypes
dSNP_nums$Rel_recessive <- dSNP_nums$nDeleterious_Hom / dSNP_nums$CalledGenotypes

#### number of heterozygous sites
dSNP_nums$Rel_NumHet <- dSNP_nums$nHet / dSNP_nums$CalledGenotypes

### proportion of heterozygous sites
dSNP_nums$PropHetGenos <- dSNP_nums$nHet / (dSNP_nums$nHet + dSNP_nums$nDeleterious_Hom)

#####################################
##### COMBINE W/ GERMPLASM INFO #####
#####################################

## merge
dSNP_nums_INFO <- merge(key[,c(9, 11, 12, 14,16,17)],
                        dSNP_nums,
                        by.x = "SequenceName",
                        by.y = "sample")

####################################
############# BOXPLOTS #############
####################################

# remove introgressed & put landraces as dots since N=3

# additive burden
boxplot1 <- Burden_boxplot(dSNP_nums_INFO,
               "group", "Rel_additive",
               "Relative Additive Burden",
               c("landrace", "OPV", "NonOil", "Oil"))

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Figure2_AdditiveBurden.pdf",
              width = 5, height = 7)
       
# smaller groupings?
Burden_boxplot(dSNP_nums_INFO,
               "Groups2", "Rel_additive",
               "Relative Additive Burden",
               c("landrace", "OPV", 
                 "HA-NonOil", "RHA-NonOil",
                 "HA-Oil", "RHA-Oil"))

# tells a similar story- reduced burden in Oil mostly seen in HA lines

# dominant burden
Burden_boxplot(dSNP_nums_INFO,
               "group", "Rel_dominant",
               "Relative Dominant Burden",
               c("landrace", "OPV", "NonOil", "Oil"))
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/DominantBurden.pdf",
       width = 5, height = 7)

# recessive burden
boxplot2 <- Burden_boxplot(dSNP_nums_INFO,
               "group", "Rel_recessive",
               "Relative Homozygous Burden",
               c("landrace", "OPV", "NonOil", "Oil"))

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Figure2_RecessiveBurden.pdf",
       width = 5, height = 7)

# number of heterozygous sites
Burden_boxplot(dSNP_nums_INFO,
               "group", "Rel_NumHet",
               "Relative Heterozygous Burden",
               c("landrace", "OPV", "NonOil", "Oil"))

# proportion of heterozygous sites
Burden_boxplot(dSNP_nums_INFO,
               "group", "PropHetGenos",
               "Proportion Heterozygous Genotypes",
               c("landrace", "OPV", "NonOil", "Oil"))

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/PropHet_AlldSNP.pdf",
       width = 5, height = 7)

####################################
###### DIFFERENCES IN MEANS? #######
####################################

# test differences in additive load between OPV, Oil, NonOil
dSNP_nums_subset <- subset(dSNP_nums_INFO, group!="landrace" &
                             group!="introgressed")
dSNP_nums_subset$group <- factor(dSNP_nums_subset$group)

# Anova to test for effects of germplasm group on additive deleterious burden
lm1 <- lm(Rel_additive ~ group,
          data=dSNP_nums_subset)
summary(lm1)
Anova(lm1) # p=0.02731
drop1(lm1, test = "F")
hist(resid(lm1))

pairs(emmeans(lm1, ~ group)) # only difference is Oil v. NonOil (p=0.0244)

# differences amoung heterotic groups?
lm2 <- lm(Rel_additive ~ group + heterotic_group +
            group:heterotic_group,
          data=dSNP_nums_subset)
Anova(lm2)
drop1(lm2, test = "F")

emmeans(lm2, ~ group|heterotic_group)
pairs(emmeans(lm2, ~ group|heterotic_group)) # diff b/w Oil and NonOil is sig for HA (p=0.0002) but not RHA (p=0.9962)

pairs(emmeans(lm2, ~ heterotic_group|group)) # diff b/w HA/RHA is sig for Oil (p=0.0001) but not NonOil (p=0.9738)


###################################### 2.) Derived dSNP data

######################################
######### SET UP DATAFRAME ###########
######################################

# file generated with Derived_Variant_Numbers.R script
SNP_genotypes <- read.table("/Users/emilydittmar/Google Drive/Active Projects/DelMutation/Results/Genotype_patterns/Annotation_VariantStats.txt", 
                            sep = "\t", header=TRUE,
                            stringsAsFactors = FALSE)

# number of derived alleles:
SNP_genotypes$NumDerivedAlleles <- 2*SNP_genotypes$NumDerivedHom + SNP_genotypes$NumHet

SNP_genotypes$Consequence <- factor(SNP_genotypes$Consequence,
                                       levels = c("SynonymousNodups", "Tolerated", "AllDel"))

######################################
############# RESHAPE ###############
######################################

# make long format:
GenoLong <- reshape(SNP_genotypes[,c(1,7,6)],
                    direction = "wide",
                    timevar = "Consequence",
                    idvar = "sample")

GenoLong$dSNP_sSNP <- GenoLong$NumDerivedAlleles.AllDel / GenoLong$NumDerivedAlleles.SynonymousNodups

Geno_all <- merge(key[,c(9, 11, 12, 14,16,17)],
                  GenoLong,
                  by.x = "SequenceName",
                  by.y = "sample")


#####################################
############# BOXPLOT ###############
#####################################

#ggplot(Geno_all, aes(x=group, y=dSNP_sSNP)) + 
#  geom_boxplot(notch = FALSE, fill="grey") +
#  theme_minimal()


# dSNP/sSNP
Burden_boxplot(Geno_all,
               "group", "dSNP_sSNP",
               "Derived dSNP/sSNP",
               c("landrace", "OPV", "NonOil", "Oil"))


ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Derived_dSNPsSNP_Boxplot.pdf",
       width = 5, height = 7)

hist(Geno_all$dSNP_sSNP)


# number of Derived deleterious alleles
Burden_boxplot(Geno_all,
               "group", "NumDerivedAlleles.AllDel",
               "Number of derived, deleterious alleles",
               c("landrace", "OPV", "NonOil", "Oil"))

######################################
######## HETERO vs. HOMOZYGOUS #######
######################################

SNP_genotypes$PropHetGenos <- SNP_genotypes$NumHet / (SNP_genotypes$NumHet+SNP_genotypes$NumDerivedHom)

#GenoProportions <- reshape(SNP_genotypes[,c(1,8,6)],
#                    direction = "wide",
#                    timevar = "Consequence",
#                    idvar = "sample")

GenoProportions <- merge(key[,c(9, 11, 12, 14, 16,17)],
                         SNP_genotypes,
                  by.x = "SequenceName",
                  by.y = "sample",
                  all.y = TRUE)

ggplot(GenoProportions, aes(x=group, y=PropHetGenos, fill = Consequence)) + 
  geom_boxplot(notch = FALSE) +
  theme_minimal()                                                                                    
## interesting that introgressed has much higher heterogyzous genotypes


boxplot3 <- ggplot(GenoProportions[which(GenoProportions$group!="landrace"),], 
       aes(x=group, y=PropHetGenos)) + 
  geom_boxplot(notch = FALSE, outlier.colour = NULL, aes(colour = Consequence)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Consequence)) +  # in order to keep outline of box black but have outlier colors the same
  theme_minimal() +
  ylab("Proportion Heterozygous Genotypes") +
  scale_x_discrete(limits=c("landrace", "OPV", "NonOil", "Oil")) +
  geom_point(data=GenoProportions[which(GenoProportions$group=="landrace" &
                                          GenoProportions$Consequence=="AllDel"),], 
             position = "identity",
             aes(x=1.25, y=PropHetGenos, fill = Consequence),
             shape=c(22,23,24), size = 3) +
  scale_fill_manual(values= c(Col_Synonymous, Col_Tolerated, Col_Deleterious),
                    name = "Variant Class", labels = c("Synonymous", "Tolerated", "Deleterious")) +
  scale_color_manual(values= c(Col_Synonymous, Col_Tolerated, Col_Deleterious),
                    name = "Variant Class", guide = "none") +
  geom_point(data=GenoProportions[which(GenoProportions$group=="landrace" &
                                          GenoProportions$Consequence=="Tolerated"),], 
             position = "identity",
             aes(x=1, y=PropHetGenos, fill = Consequence),
             shape=c(22,23,24), size = 3) +
  geom_point(data=GenoProportions[which(GenoProportions$group=="landrace" &
                                          GenoProportions$Consequence=="SynonymousNodups"),], 
             position = "identity",
             aes(x=0.75, y=PropHetGenos, fill = Consequence),
             shape=c(22,23,24), size = 3) +
  theme(axis.text = element_text(size=10),
        axis.text.x = element_text(angle=90),
        #legend.position = "none",
        legend.text = element_text(size=12),
        strip.text.x = element_text(size = 12))

# save with and without legend to enlarge  
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Figure2_GenotypeBoxplot_noLeg.pdf",
       width = 5, height = 7)
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Figure2_GenotypeBoxplot_wLeg.pdf",
       width = 5, height = 7)

######################################
####### COMBINE 3 FOR MAIN TEXT ######
######################################

ggarrange(boxplot1, boxplot2, boxplot3, nrow=1)
