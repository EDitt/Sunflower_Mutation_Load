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

key <- read.csv("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Sunflower_Mutation_Load/BAD_Mutations/LineKeywINFO.csv", 
                header=T)

######################################
############# RESHAPE ###############
######################################


# make long format:
GenoLong <- reshape(SNP_genotypes[,c(1,7,6)],
                    direction = "wide",
                    timevar = "Consequence",
                    idvar = "sample")

GenoLong$dSNP_sSNP <- GenoLong$NumDerivedAlleles.AllDel / GenoLong$NumDerivedAlleles.SynonymousNodups

Geno_all <- merge(key[c(9, 11, 12, 14)],
                  GenoLong,
                  by.x = "SequenceName",
                  by.y = "sample")

#####################################
###### ADJUST CATEGORIZATIONS #######
#####################################

# separate landrace & OPV
# combine HA/RHA- keep as Oil vs. NonOil

Geno_all$group <- ifelse(Geno_all$heterotic_group=="landrace",
                         "landrace",
                         ifelse(Geno_all$heterotic_group=="HA" |
                                  Geno_all$heterotic_group=="RHA",
                                as.character(Geno_all$Oil_NonOil),
                                as.character(Geno_all$Groups)))

Geno_all$group[which(Geno_all$group=="landrace_OPV")] <- "OPV"

aggregate(Geno_all$SequenceName, by=list(Geno_all$group, Geno_all$Groups), length) # check


levels(as.factor(Geno_all$group))
Geno_all$group <- factor(Geno_all$group, levels=c("landrace", "OPV", "introgressed",
                                                    "NonOil", "Oil"))

#Geno_all$group <- factor(Geno_all$group, levels=c("landrace", "Oil", "introgressed",
#                                                  "OPV", "NonOil"))

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


