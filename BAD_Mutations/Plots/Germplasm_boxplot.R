## derived dSNPs/sSNPs for groups of germplasm

source("BAD_Mutations/Variant_analyses/Functions.R")

# formerly Genotype_Load_Analyses.R

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

key <- read.csv("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Sunflower_Mutation_Load/BAD_Mutations/LineKeywINFO.csv", 
                header=T)

#####################################
###### ADJUST CATEGORIZATIONS #######
#####################################

# separate landrace & OPV
# combine HA/RHA- keep as Oil vs. NonOil

key$group <- ifelse(key$heterotic_group=="landrace",
                         "landrace",
                         ifelse(key$heterotic_group=="HA" |
                                  key$heterotic_group=="RHA",
                                as.character(key$Oil_NonOil),
                                as.character(key$Groups)))

key$group[which(key$group=="landrace_OPV")] <- "OPV"


aggregate(key$sort, by=list(key$group, key$Groups), length) # check


levels(as.factor(key$group))
key$group <- factor(key$group, levels=c("landrace", "OPV", "introgressed",
                                                  "NonOil", "Oil"))

######################################
############# RESHAPE ###############
######################################


# make long format:
GenoLong <- reshape(SNP_genotypes[,c(1,7,6)],
                    direction = "wide",
                    timevar = "Consequence",
                    idvar = "sample")

GenoLong$dSNP_sSNP <- GenoLong$NumDerivedAlleles.AllDel / GenoLong$NumDerivedAlleles.SynonymousNodups

Geno_all <- merge(key[,c(9, 11, 12, 14,16)],
                  GenoLong,
                  by.x = "SequenceName",
                  by.y = "sample")


#####################################
############# BOXPLOT ###############
#####################################

ggplot(Geno_all, aes(x=group, y=dSNP_sSNP)) + 
  geom_boxplot(notch = FALSE, fill="grey") +
  theme_minimal()

# remove introgressed & put landraces as dots since N=3

ggplot(Geno_all[which(Geno_all$group!="landrace"),], aes(x=group, y=dSNP_sSNP)) + 
  geom_boxplot(notch = FALSE, fill="grey") +
  theme_minimal() +
  ylab("derived dSNPs/sSNPs") +
  scale_x_discrete(limits=c("landrace", "OPV", "NonOil", "Oil")) +
  geom_point(data=Geno_all[which(Geno_all$group=="landrace"),], aes(y=dSNP_sSNP),
             shape=8, size = 2) +
  theme(axis.text = element_text(size=10),
        axis.text.x = element_text(angle=90),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size = 12))

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/Figure2_Boxplot.pdf",
       width = 5, height = 7)

hist(Geno_all$dSNP_sSNP)


######################################
######## HETERO vs. HOMOZYGOUS #######
######################################

SNP_genotypes$PropHetGenos <- SNP_genotypes$NumHet / (SNP_genotypes$NumHet+SNP_genotypes$NumDerivedHom)

#GenoProportions <- reshape(SNP_genotypes[,c(1,8,6)],
#                    direction = "wide",
#                    timevar = "Consequence",
#                    idvar = "sample")

GenoProportions <- merge(key[,c(9, 11, 12, 14, 16)],
                         SNP_genotypes,
                  by.x = "SequenceName",
                  by.y = "sample",
                  all.y = TRUE)

ggplot(GenoProportions, aes(x=group, y=PropHetGenos, fill = Consequence)) + 
  geom_boxplot(notch = FALSE) +
  theme_minimal()                                                                                    
## interesting that introgressed has much higher heterogyzous genotypes


ggplot(GenoProportions[which(GenoProportions$group!="landrace"),], 
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
             shape=c(21,22,23), size = 3) +
  scale_fill_manual(values= c(Col_Synonymous, Col_Tolerated, Col_Deleterious),
                    name = "Variant Class", labels = c("Synonymous", "Tolerated", "Deleterious")) +
  scale_color_manual(values= c(Col_Synonymous, Col_Tolerated, Col_Deleterious),
                    name = "Variant Class", guide = "none") +
  geom_point(data=GenoProportions[which(GenoProportions$group=="landrace" &
                                          GenoProportions$Consequence=="Tolerated"),], 
             position = "identity",
             aes(x=1, y=PropHetGenos, fill = Consequence),
             shape=c(21,22,23), size = 3) +
  geom_point(data=GenoProportions[which(GenoProportions$group=="landrace" &
                                          GenoProportions$Consequence=="SynonymousNodups"),], 
             position = "identity",
             aes(x=0.75, y=PropHetGenos, fill = Consequence),
             shape=c(21,22,23), size = 3) +
  theme(axis.text = element_text(size=10),
        axis.text.x = element_text(angle=90),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size = 12))
  


