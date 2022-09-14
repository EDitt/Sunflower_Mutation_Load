### Comparing inside/outside hotspots of LROH

source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")
source("BAD_Mutations/Variant_analyses/Functions.R")

### 1.) Frequency distribution inside/outside of ROH hotspots

######################################
######### SET UP DATAFRAME ###########
######################################

IN_sfs <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/DerivedFreq_INhotspots_Bins.txt",
                           sep="\t", header = T)
IN_sfs$roh <- "insideROH"
OUT_sfs <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/DerivedFreq_OUThotspots_Bins.txt",
                     sep="\t", header = T)
OUT_sfs$roh <- "outsideROH"

all_sfs <- rbind(IN_sfs, OUT_sfs)

AnnotationsToKeep <- c("SynonymousNodups", "AllDel")

sfs_subset <- subset(all_sfs, Annotation %in% AnnotationsToKeep)
sfs_subset$Annotation <- factor(sfs_subset$Annotation, levels=AnnotationsToKeep)
sfs_subset$roh <- factor(sfs_subset$roh, levels=c("insideROH", "outsideROH"))

ggplot(sfs_subset, aes(x=breaks, y=prop, fill=Annotation)) +
  geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  scale_fill_manual(values= c(Col_Synonymous, Col_Deleterious),
                    name = "Variant Class", labels = c("Synonymous", "Deleterious")) +
  ylab("Proportion") +
  xlab("Frequency") +
  facet_wrap(~roh) +
  theme(panel.grid.minor = element_blank())

ggplot(sfs_subset, aes(x=breaks, y=prop, fill=roh)) +
  geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  #scale_fill_manual(values= c(Col_Synonymous, Col_Deleterious),
  #                  name = "Variant Class", labels = c("Synonymous", "Deleterious")) +
  ylab("Proportion") +
  xlab("Frequency") +
  facet_wrap(~Annotation) +
  theme(panel.grid.minor = element_blank())

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/SFS_ROH.pdf")

  #scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
  #                   labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", 
  #                              "[0.20,0.25)", "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)")) +

######################################
##### DENSITY PLOT? (ON CLUSTER) ######
######################################

IN_ROH <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Hotspots/SNPs_in_rohHotspots_HEADER.txt",
                     header=T)
IN_ROH$ROH <- "ROH"

OUT_ROH <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Hotspots/SNPs_outside_rohHotspots_HEADER.txt",
                      header=T)
OUT_ROH$ROH <- "OUTSIDE_ROH"

Both_ROH <- rbind(IN_ROH, OUT_ROH)

Both_ROH$ROH <- factor(Both_ROH$ROH, levels=c("ROH", "OUTSIDE_ROH"))

#Both_ROH_del <- subset(Both_ROH, Variant_type=="AllDel")

plot <- ggplot(Both_ROH[which(Both_ROH$Variant_type=="AllDel"),], aes(x=Derived_Freq)) + 
  geom_density(alpha=.6, aes(fill=ROH)) +
  theme_minimal() +
  xlab("Derived Frequency")

pdf(file="/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Hotspots/dSNP_density.pdf")
print(plot)
dev.off()

plot2 <- ggplot(Both_ROH[which(Both_ROH$Variant_type=="SynonymousNodups"),], aes(x=Derived_Freq)) + 
  geom_density(alpha=.6, aes(fill=ROH)) +
  theme_minimal() +
  xlab("Derived Frequency")

pdf(file="/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Hotspots/sSNP_density.pdf")
print(plot2)
dev.off()

### 2.) Number of dSNPs/codon, sSNPs/codon, dSNPs/sSNPs

######################################
######### SET UP DATAFRAME ###########
######################################

# combine numbers of codons, dSNPs, sSNPs in/out of LROH

ROH_info <- ImportFilesAsList("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH",
                              "_LROH.bed", "Num_", FALSE)

for (i in seq_along(ROH_info)) {
  name <- names(ROH_info[i])
  colnames(ROH_info[[i]]) <- c("Chromosome", "StartPos", "EndPos", paste0("NumVariants_", name))
}

Codon_info <- ImportFilesAsList("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH",
                                "_NumCodonCounts.txt", "", FALSE)

for (i in seq_along(Codon_info)) {
  name <- names(Codon_info[i])
  colnames(Codon_info[[i]]) <- c("Chromosome", "StartPos", "EndPos", "TotalCodingBases", "TotalmRNA", "Num_Codons")
}

# In ROH: (N=316)
ROH_data <- merge(Codon_info$ROH,
                  merge(ROH_info$AllDel_IN, ROH_info$SynonymousNodups_IN, by=c("Chromosome", "StartPos", "EndPos")),
                  by=c("Chromosome", "StartPos", "EndPos"), all=TRUE)

ROH_data$dSNP_codon <- ROH_data$NumVariants_AllDel_IN / ROH_data$Num_Codons
ROH_data$sSNP_codon <- ROH_data$NumVariants_SynonymousNodups_IN / ROH_data$Num_Codons
ROH_data$dSNP_sSNP <- ROH_data$NumVariants_AllDel_IN / ROH_data$NumVariants_SynonymousNodups_IN
ROH_data$ROH_Cat <- "In"

ROH_data_noNA <- subset(ROH_data, NumVariants_AllDel_IN > 0 | NumVariants_SynonymousNodups_IN > 0 | !is.na(Num_Codons))
length(ROH_data_noNA$Chromosome) #255

# Outside ROH: # N=12409
OUT_ROH_data <- merge(Codon_info$OUT_ROH,
                      merge(ROH_info$AllDel_OUT, ROH_info$SynonymousNodups_OUT, by=c("Chromosome", "StartPos", "EndPos")),
                      by=c("Chromosome", "StartPos", "EndPos"), all.x=TRUE, all.y=TRUE)

OUT_ROH_data$dSNP_codon <- OUT_ROH_data$NumVariants_AllDel_OUT / OUT_ROH_data$Num_Codons
OUT_ROH_data$sSNP_codon <- OUT_ROH_data$NumVariants_SynonymousNodups_OUT / OUT_ROH_data$Num_Codons
OUT_ROH_data$dSNP_sSNP <- OUT_ROH_data$NumVariants_AllDel_OUT / OUT_ROH_data$NumVariants_SynonymousNodups_OUT
OUT_ROH_data$ROH_Cat <- "Out"

OUT_ROH_data_noNA <- subset(OUT_ROH_data, NumVariants_AllDel_OUT > 0 | NumVariants_SynonymousNodups_OUT > 0 | !is.na(Num_Codons))
length(OUT_ROH_data_noNA$Chromosome) # 11300

All_ROH_data <- rbind(ROH_data_noNA [,c("Chromosome", "StartPos", "EndPos", "dSNP_codon", "sSNP_codon", "dSNP_sSNP", "ROH_Cat")],
                      OUT_ROH_data_noNA[,c("Chromosome", "StartPos", "EndPos", "dSNP_codon", "sSNP_codon", "dSNP_sSNP", "ROH_Cat")])

### saved on cluster
#write.table(All_ROH_data, file="/scratch/eld72413/SAM_seq/dSNP_results/GenomicPatterns/LROH/Hotspots/ROH_VariantNums.txt",
#            row.names=FALSE, quote=FALSE, sep="\t")
### load locally
#ROH_VarNums <- read.table(file="/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/ROH/ROH_VariantNums.txt",
#                          header=T)



# change NA's to 0 (sSNP num is zero)

All_ROH_data$dSNP_sSNP_noNA <- ifelse(is.na(All_ROH_data$dSNP_sSNP) | is.infinite(All_ROH_data$dSNP_sSNP),
                                     0, All_ROH_data$dSNP_sSNP)

All_ROH_data[which(is.na(All_ROH_data$dSNP_codon)),]
All_ROH_data[which(is.na(All_ROH_data$sSNP_codon)),]

#ROH_VarNums$dSNP_sSNP[which(is.na(ROH_VarNums$dSNP_sSNP))] <- 0
#ROH_VarNums$dSNP_sSNP[which(is.infinite(ROH_VarNums$dSNP_sSNP))] <- 0
hist(All_ROH_data$dSNP_sSNP)
hist(All_ROH_data$dSNP_sSNP_noNA)


######################################
########### Num Per Codon ############
######################################

ROH_VarNums_long <- reshape(All_ROH_data,
                            varying=list(names(ROH_VarNums)[c(4:5)]),
                            v.names = "SNP_codon",
                            idvar=c("Chromosome", "StartPos", "EndPos", "ROH_Cat"),
                            timevar = "Variant_type",
                            times = names(ROH_VarNums)[c(4:5)],
                            drop=names(ROH_VarNums)[c(6,8)],
                            direction="long"
)

ggplot(data=ROH_VarNums_long, aes(x=ROH_Cat, y=SNP_codon)) +
  geom_boxplot(aes(color=Variant_type))

ggplot(data=ROH_VarNums_long, aes(x=Variant_type, y=SNP_codon)) +
  geom_boxplot()+
  facet_wrap(~ROH_Cat, scales = "free_y") # a bit misleading- this one is better for showing dSNP/sSNP ratio patterns

### ***** this is probably the best way to portray:
ggplot(data=ROH_VarNums_long, aes(x=ROH_Cat, y=SNP_codon)) +
  geom_boxplot()+
  #geom_violin() +
  facet_wrap(~Variant_type, scales = "free_y")

min(ROH_VarNums_long[which(ROH_VarNums_long$SNP_codon!=0), "SNP_codon"]) # 0.00003245
# add constant so zero values aren't removed?

ggplot(data=ROH_VarNums_long, aes(x=ROH_Cat, y=log(SNP_codon))) +
  geom_boxplot()+
  facet_wrap(~Variant_type, scales = "free_y")

ggplot(data=ROH_VarNums_long, aes(x=ROH_Cat, y=log(SNP_codon + 0.000001))) +
  geom_boxplot()+
  #geom_violin() +
  facet_wrap(~Variant_type, scales = "free_y")

ggplot(ROH_VarNums, aes(x=dSNP_codon)) + 
  geom_density(alpha=.6, aes(fill=ROH_Cat)) +
  theme_minimal() +
  xlim(0,0.005) +
  xlab("dSNP/codon")
ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/ROH_dSNPcodon.pdf")

ROH_VarNums[which(ROH_VarNums$dSNP_codon>0.005 &
                    ROH_VarNums$ROH_Cat=="In"),] # Chrom 5,12,17

######################################
############# dSNP/sSNP ##############
######################################

ROH_VarNums <- All_ROH_data

min(ROH_VarNums[which(ROH_VarNums$dSNP_sSNP!=0), "dSNP_sSNP"]) 

ggplot(data=ROH_VarNums, aes(x=ROH_Cat, y=log(dSNP_sSNP + 0.0001))) +
  geom_boxplot()
  #geom_violin()

ggplot(data=ROH_VarNums, aes(x=ROH_Cat, y=log(dSNP_sSNP))) +
  geom_boxplot() # removing zeros keeps them relatively equal

ggplot(data=ROH_VarNums, aes(x=ROH_Cat, y=dSNP_sSNP)) +
  geom_boxplot()
  #geom_violin()

ggplot(data=ROH_VarNums, aes(dSNP_sSNP, fill=ROH_Cat)) +
  geom_histogram(binwidth = 0.1, alpha=0.8)

ggplot(ROH_VarNums, aes(x=dSNP_sSNP)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.1) +
  #geom_density(alpha=.2, fill="#FF6666") +
  geom_density(alpha=.6, aes(fill=ROH_Cat)) +
  theme_minimal()

# ****this is probably the best way to see distribution-
ggplot(ROH_VarNums, aes(x=dSNP_sSNP_noNA)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.1) +
  #geom_density(alpha=.2, fill="#FF6666") +
  geom_density(alpha=.6, aes(fill=ROH_Cat)) +
  theme_minimal() +
  xlim(0,0.75) +
  xlab("dSNP/sSNP")

ggsave("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Manuscript/Sunflower_MutationLoad_Manuscript/Sunflower_MutationLoad_v1/RawFigs/ROH_dSNPsSNP.pdf")

ggplot(ROH_VarNums, aes(x=log(dSNP_sSNP))) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.1) +
  #geom_density(alpha=.2, fill="#FF6666") +
  geom_density(alpha=.6, aes(fill=ROH_Cat)) +
  theme_minimal()

ggplot(ROH_VarNums, aes(x=log(dSNP_sSNP + 0.0001))) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.1) +
  #geom_density(alpha=.2, fill="#FF6666") +
  geom_density(alpha=.6, aes(fill=ROH_Cat)) +
  theme_minimal()

ggplot(ROH_VarNums, aes(x=log(dSNP_sSNP_noNA + 0.0001))) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.1) +
  #geom_density(alpha=.2, fill="#FF6666") +
  geom_density(alpha=.6, aes(fill=ROH_Cat)) +
  theme_minimal()

ROH_VarNums$Length <- ROH_VarNums$EndPos - ROH_VarNums$StartPos
hist(ROH_VarNums$Length)
max(ROH_VarNums$Length)

ggplot(ROH_VarNums, aes(x=Length, y=dSNP_sSNP_noNA)) +
  geom_smooth(method="lm", aes(color=ROH_Cat))

ggplot(ROH_VarNums, aes(x=Length, y=dSNP_sSNP)) +
  geom_smooth(method="lm", aes(color=ROH_Cat)) +
  geom_point(aes(color=ROH_Cat))

ggplot(ROH_VarNums, aes(x=Length, y=dSNP_sSNP)) +
  geom_smooth(method="lm") +
  geom_point() +
  facet_wrap(~ROH_Cat, scales = "free")


