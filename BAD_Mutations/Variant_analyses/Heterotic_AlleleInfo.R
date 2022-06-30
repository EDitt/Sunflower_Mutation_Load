# Data Wrangling for allelic differentiation between heterotic groups
### 1.) alleles that are private to HA/RHA
### 2.) frequency differentiation

source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/Functions.R")


######################################
######### SET UP DATAFRAMES ##########
######################################

Groups <- ImportFilesAsList("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup", ".txt", "_SNP_info_ForUnfolded_reduced")
Groups <- lapply(Groups, function(x) {
	x$Num_Ref_alleles <- x$Num_alleles - x$Num_Alt_alleles; return(x)
	}
	)

Oil <- merge(Groups[["HA-Oil"]], Groups[["RHA-Oil"]], 
	by=c("Chromosome", "Position", "Ref_allele", "Alt_allele", "Alt_Freq", "Ancestral_Allele", "Variant_type"),
	suffixes=c("_HA", "_RHA"))
NonOil <- merge(Groups[["HA-NonOil"]], Groups[["RHA-NonOil"]], 
	by=c("Chromosome", "Position", "Ref_allele", "Alt_allele", "Alt_Freq", "Ancestral_Allele", "Variant_type"),
	suffixes=c("_HA", "_RHA"))

All <- merge(Oil, NonOil, 
	by=c("Chromosome", "Position", "Ref_allele", "Alt_allele", "Alt_Freq", "Ancestral_Allele", "Variant_type"),
	suffixes=c("Oil", "NonOil"))

All$MAF <- ifelse(All$Alt_Freq < 0.05, All$Alt_Freq,
	1-All$Alt_Freq)

# check
length(All[which(All$Alt_Freq!=All$MAF),"Position"]) # 344,050
length(All[which(All$Alt_Freq==All$MAF),"Position"]) # 665,280

# remove alleles that are invariant across all lines
All_subset <- subset(All, Num_Alt_alleles_HAOil + Num_Alt_alleles_RHAOil + Num_Alt_alleles_HANonOil + Num_Alt_alleles_RHANonOil > 0 &
	Num_Ref_alleles_HAOil + Num_Ref_alleles_RHAOil + Num_Ref_alleles_HANonOil + Num_Ref_alleles_RHANonOil > 0 
	) # N= 825,503 (825,869 when removing invariant alternate allele positions)
length(All_subset[which(All_subset$Alt_Freq!=All_subset$MAF),"Position"]) # 343,684


#####################################
##### MERGE INFO - ALL ALLELES ######
#####################################

# combine Oil and NonOil
All_subset$Num_Alt_alleles_HA <- All_subset$Num_Alt_alleles_HAOil + All_subset$Num_Alt_alleles_HANonOil
All_subset$Num_Ref_alleles_HA <- All_subset$Num_Ref_alleles_HAOil + All_subset$Num_Ref_alleles_HANonOil
All_subset$Num_Alt_alleles_RHA <- All_subset$Num_Alt_alleles_RHAOil + All_subset$Num_Alt_alleles_RHANonOil
All_subset$Num_Ref_alleles_RHA <- All_subset$Num_Ref_alleles_RHAOil + All_subset$Num_Ref_alleles_RHANonOil

write.table(All_subset, "/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo.txt", 
	sep = "\t", quote=FALSE, row.names=FALSE)

# HA v. RHA for Oil, NonOil, across both groups
Oil_FreqBins <- Group_freqbins("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo.txt", 
	"HAOil", "RHAOil", "Num_Alt_alleles_", "Num_Ref_alleles_", 0.5, 0.1, "MAF")
NonOil_FreqBins <- Group_freqbins("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo.txt", 
	"HANonOil", "RHANonOil", "Num_Alt_alleles_", "Num_Ref_alleles_", 0.5, 0.1, "MAF")
Both_FreqBins <- Group_freqbins("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo.txt", 
	"HA", "RHA", "Num_Alt_alleles_", "Num_Ref_alleles_", 0.5, 0.1, "MAF")

#####################################
####### MERGE INFO - DERIVED ########
#####################################

# look only at derived
All_derived <- subset(All_subset, Ancestral_Allele!="NA") # N= 593,187
All_derived$Overall_Derived_Freq <- ifelse(All_derived$Ref_allele==All_derived$Ancestral_Allele,
	All_derived$Alt_Freq,
	ifelse(All_derived$Alt_allele==All_derived$Ancestral_Allele,
		1-All_derived$Alt_Freq,
		NA))
length(All_derived[which(is.na(All_derived$Overall_Derived_Freq)),]) # 41- both derived

write.table(All_derived, "/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo_derived.txt", 
	sep = "\t", quote=FALSE, row.names=FALSE)

# HA v. RHA for Oil, NonOil, across both groups - derived
Oil_DerivedFreqBins <- Group_freqbins("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo_derived.txt", 
	"HAOil", "RHAOil", "Num_Alt_alleles_", "Num_Ref_alleles_", 1.0, 0.1, "Overall_Derived_Freq")
NonOil_DerivedFreqBins <- Group_freqbins("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo_derived.txt", 
	"HANonOil", "RHANonOil", "Num_Alt_alleles_", "Num_Ref_alleles_", 1.0, 0.1, "Overall_Derived_Freq")
Both_DerivedFreqBins <- Group_freqbins("/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo_derived.txt", 
	"HA", "RHA", "Num_Alt_alleles_", "Num_Ref_alleles_", 1.0, 0.1, "Overall_Derived_Freq")
# 588 private alleles in the 0.9-1.0 derived frequency class

HighFreqDerived <- subset(All_derived,
	Variant_type=="AllDel" &
	Overall_Derived_Freq >= 0.9 &
	((Num_Alt_alleles_RHA > 1 & Num_Alt_alleles_HA==0) |
	(Num_Alt_alleles_RHA ==0 & Num_Alt_alleles_HA>1))) # 584


Chromosomes <- aggregate(HighFreqDerived$Position, by=list(HighFreqDerived$Chromosome), length) # spread evenly?

# save to local computer
save(Oil_FreqBins, NonOil_FreqBins, Both_FreqBins, 
	Oil_DerivedFreqBins, NonOil_DerivedFreqBins, Both_DerivedFreqBins,
	file="FreqBins.RData")

# also using "/scratch/eld72413/SAM_seq/dSNP_results/GenotypeInfo/HeteroticGroup/AlleleInfo_derived.txt" for joint SFS figures