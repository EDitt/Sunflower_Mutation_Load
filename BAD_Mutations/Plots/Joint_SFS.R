### code for Joint SFS figures

source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")

# plot code formerly in README.md

######################################
######### SET UP DATAFRAME ###########
######################################


snp_freqs_derived <- read.table("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo_derived.txt", header=T, sep="\t")


snp_freqs_type <- split(snp_freqs_derived, snp_freqs_derived$Variant_type)


######################################
############ OIL LINES ##############
######################################

SFS_Oilplots1 <- lapply(snp_freqs_type, function(x) {
	JointSFSPlot1(x, 0.05, "Derived_Freq_HAOil", "Derived_Freq_RHAOil", "HA-Oil Derived Freq", "RHA-Oil Derived Freq")
	})
pdf("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/HA_RHA_OilJointSFS.pdf")
print(ggarrange(SFS_Oilplots1$SynonymousNodups, SFS_Oilplots1$AllDel), labels=c("A. Synonymous", "B. Deleterious"))
dev.off()

SFS_Oilplots2 <- lapply(snp_freqs_type, function(x) {
	JointSFSPlot2(x, 100, "Derived_Freq_HAOil", "Derived_Freq_RHAOil", "HA-Oil Derived Freq", "RHA-Oil Derived Freq")
	})
pdf("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/HA_RHA_OilJointSFS2.pdf")
print(ggarrange(SFS_Oilplots2$SynonymousNodups, SFS_Oilplots2$AllDel), labels=c("A. Synonymous", "B. Deleterious"))
dev.off()

######################################
########### NONOIL LINES #############
######################################

SFS_NonOilplots1 <- lapply(snp_freqs_type, function(x) {
	JointSFSPlot1(x, 0.05, "Derived_Freq_HANonOil", "Derived_Freq_RHANonOil", "HA-NonOil Derived Freq", "RHA-NonOil Derived Freq")
	})
pdf("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/HA_RHA_NonOilJointSFS.pdf")
print(ggarrange(SFS_NonOilplots1$SynonymousNodups, SFS_NonOilplots1$AllDel), labels=c("A. Synonymous", "B. Deleterious"))
dev.off()

SFS_NonOilplots2 <- lapply(snp_freqs_type, function(x) {
	JointSFSPlot2(x, 100, "Derived_Freq_HANonOil", "Derived_Freq_RHANonOil", "HA-NonOil Derived Freq", "RHA-NonOil Derived Freq")
	})
pdf("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/HA_RHA_NonOilJointSFS2.pdf")
print(ggarrange(SFS_NonOilplots2$SynonymousNodups, SFS_NonOilplots2$AllDel), labels=c("A. Synonymous", "B. Deleterious"))
dev.off()

######################################
############# COMBINED ###############
######################################

# wrangling to get overall HA and RHA frequencies:
snp_freqs_derived$Derived_Freq_HA <- ifelse(snp_freqs_derived$Ref_allele==snp_freqs_derived$Ancestral_Allele,
	snp_freqs_derived$Num_Alt_alleles_HA / (snp_freqs_derived$Num_Alt_alleles_HA + snp_freqs_derived$Num_Ref_alleles_HA),
	ifelse(snp_freqs_derived$Alt_allele==snp_freqs_derived$Ancestral_Allele,
		snp_freqs_derived$Num_Ref_alleles_HA / (snp_freqs_derived$Num_Alt_alleles_HA + snp_freqs_derived$Num_Ref_alleles_HA),
		NA))

snp_freqs_derived$Derived_Freq_RHA <- ifelse(snp_freqs_derived$Ref_allele==snp_freqs_derived$Ancestral_Allele,
	snp_freqs_derived$Num_Alt_alleles_RHA / (snp_freqs_derived$Num_Alt_alleles_RHA + snp_freqs_derived$Num_Ref_alleles_RHA),
	ifelse(snp_freqs_derived$Alt_allele==snp_freqs_derived$Ancestral_Allele,
		snp_freqs_derived$Num_Ref_alleles_RHA / (snp_freqs_derived$Num_Alt_alleles_RHA + snp_freqs_derived$Num_Ref_alleles_RHA),
		NA))

snp_freqs_type <- split(snp_freqs_derived, snp_freqs_derived$Variant_type) # re-run split command

SFS_plots1 <- lapply(snp_freqs_type, function(x) {
	JointSFSPlot1(x, 0.05, "Derived_Freq_HA", "Derived_Freq_RHA", "HA Derived Freq", "RHA Derived Freq")
	})
pdf("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/HA_RHA_JointSFS.pdf")
print(ggarrange(SFS_plots1$SynonymousNodups, SFS_plots1$AllDel), labels=c("A. Synonymous", "B. Deleterious"))
dev.off()

SFS_plots2 <- lapply(snp_freqs_type, function(x) {
	JointSFSPlot2(x, 100, "Derived_Freq_HA", "Derived_Freq_RHA", "HA Derived Freq", "RHA Derived Freq")
	})
pdf("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/HA_RHA_JointSFS2.pdf")
print(ggarrange(SFS_plots2$SynonymousNodups, SFS_plots2$AllDel), labels=c("A. Synonymous", "B. Deleterious"))
dev.off()




