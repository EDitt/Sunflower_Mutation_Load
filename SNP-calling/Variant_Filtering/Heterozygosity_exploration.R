setwd("/Users/eld72413/Google Drive/Active Projects/DelMutation/Variant_graphs")

hets_full <- read.csv("Variants_HetInfo.table", sep = "\t", header = TRUE)
length(hets_full$CHROM) #4,903,270
hets_full$het_prop <- hets_full$HET / hets_full$NCALLED

hist(hets_full$het_prop)

# account for NCALLED > 230
hets <- subset(hets_full, NCALLED >230)
length(hets$CHROM) #3,848,922

########### how many sites have a greater than 20% proportion of heterozygotes?
length(which(hets$het_prop > 0.2)) # 142,975 (3.7%)
# 166,638 (3.4%) out of full set
length(which(hets$HET > 57)) #134,714
# 144,314 out of full set

##### why didn't these get filtered out?
head(hets[which(hets$het_prop > 0.2),])
### I used grep to see if there was any difference in these sites before and after the 'SelectVariants' script. There was not

##### Use ExcessHet filter?
hist(hets$ExcessHet)
hist(hets[which(hets$ExcessHet < 20),"ExcessHet"])
hist(hets[which(hets$ExcessHet < 10),"ExcessHet"])

### what is the overlap between high excess het and high het. proportion?
HighHetProp <- which(hets$het_prop > 0.2)
HighExcessHet <- which(hets$ExcessHet > 5)
length(HighExcessHet) # 382,168
#395,838 (8%) when looking at full set
length(which(HighHetProp %in% HighExcessHet)) # 108,118 (76% of het prop > 0.2)
#112,639 when looking across full set

### which ones are not included? 
HighHetPropLowExcessHet <-hets[which(hets$het_prop > 0.2 & hets$ExcessHet < 5),]
length(HighHetPropLowExcessHet$CHROM) #34,857
head(HighHetPropLowExcessHet) #high numbers of homozygous variant
hist(HighHetPropLowExcessHet$HET)

#Filter by ExcessHet greater than 5?
ExcessHet5_filter <- subset(hets, ExcessHet < 5)
# smallest number of heterozygotes?
max(ExcessHet5_filter$HET) #285
hist(ExcessHet5_filter$ExcessHet)
hist(ExcessHet5_filter$HET)
head(ExcessHet5_filter[which(ExcessHet5_filter$HET == 285),]) # Excess Het of 4.3474

# what did I remove?
HighExcessHet5 <- subset(hets, ExcessHet > 5)
min(HighExcessHet5$HET)
head(HighExcessHet5[which(HighExcessHet5$HET == 0),])
length(HighExcessHet5[which(HighExcessHet5$HET == 0),"CHROM"]) #884 have 0 hets???

hist(hets$InbreedingCoeff)


########## AFTER USING BCFTOOLS FILTERS
# test to make sure it filtered out sites with more than 20% heterozygotes
hets2 <- read.csv("bcftools_Variants_HetInfo.table", sep = "\t", header = TRUE)

hets2$het_prop <- hets2$HET / hets2$NCALLED
hist(hets2$het_prop)
max(hets2$het_prop) #0.1993007

hist(hets2$ExcessHet)
hist(hets2[which(hets2$ExcessHet < 20),"ExcessHet"])
hist(hets2[which(hets2$ExcessHet < 10),"ExcessHet"])

hets2_sub <- subset(hets2, NCALLED > 230) #N = 3,705,353
head(hets2_sub[which(hets2_sub$ExcessHet > 10),])
head(hets2_sub[which(hets2_sub$ExcessHet > 5),])
length(hets2_sub[which(hets2_sub$ExcessHet > 10),"CHROM"]) #100,890
length(hets2_sub[which(hets2_sub$ExcessHet > 5),"CHROM"]) #273,693

## what is the highest ExcessHet value for sites with only 1 heterozygote?
hets2_1Het <- subset(hets2_sub, HET==1)
max(hets2_1Het$ExcessHet) #91.56
head(hets2_1Het[which(hets2_1Het$ExcessHet > 75),]) #N=12
hist(hets2_1Het$ExcessHet)

# filter at ExcessHet=5?
hets2_5EH <- subset(hets2_sub, ExcessHet < 5)
length(hets2_5EH$CHROM) #3,431,660 out of 3,705,353 (92.6%)
hist(hets2_5EH$ExcessHet)
hist(hets2_5EH$InbreedingCoeff)
hist(hets2_5EH$het_prop)
hist(hets2_sub$het_prop)

# what did I filter out
hets2_5EHfilt <- subset(hets2_sub, ExcessHet > 5)
hist(hets2_5EHfilt$het_prop)
hist(hets2_5EH$het_prop)
hist(hets2_sub$het_prop)

hist(hets2_5EH$InbreedingCoeff)
hist(hets2_5EHfilt$InbreedingCoeff)

# how many have ExcessHet values greater than 5? 10?
length(hets2_1Het[which(hets2_1Het$ExcessHet > 5),"CHROM"]) #7842
length(hets2_1Het[which(hets2_1Het$ExcessHet > 10),"CHROM"]) #2381

### what about using the inbreeding coefficient? GATK says high values can be a proxy for poor mapping
# positive numbers indicate fewer than expected heterozygotes which imply inbreeding (which these samples are!)
hist(hets2_sub$InbreedingCoeff)
head(hets2_sub[which(hets2_sub$InbreedingCoef>0.6),])

head(hets2_sub[which(hets2_sub$InbreedingCoef<0),])

head(hets2_sub[which(hets2_sub$InbreedingCoef < -0.3),])

# smallest number for 1 het?
min(hets2_1Het$InbreedingCoeff) # -0.265

min(na.omit(hets2_sub$InbreedingCoeff)) #-0.393

### filter out sites with Inbreeding coefficient less than -0.3?
hets2_filterIC <- subset(hets2_sub, InbreedingCoeff < -0.3)
length(hets2_filterIC$CHROM) #14
hist(hets2_filterIC$ExcessHet)

######################### SCRATCH
head(hets[HighExcessHet,])
hist(hets[HighExcessHet,"het_prop"])
length(which(hets$ExcessHet > 5 & hets$het_prop < 0.2)) #282,824
head(hets[which(hets$ExcessHet > 5 & hets$het_prop < 0.2),])


# Distribution of ExcessHet among samples with > 20% heterozygous
hets_highHetProp <- subset(hets, het_prop > 0.2)
hist(hets_highHetProp$ExcessHet)
hist(hets$ExcessHet)

hets_highHetProp2 <- subset(hets, het_prop > 0.3)
hist(hets_highHetProp2$ExcessHet) 

# "ExcessHet" takes into account the # of homozygous non-reference genotypes
# A locus with 0.24 % heterozygous samples has 54/229 homozygous non-reference

### distribution of heterzygote proportion with high excess het
HighExcessHet2 <- which(hets$ExcessHet > 2)
length(HighExcessHet2) #2,024,868 (41%)

hets_highexcesshet <- subset(hets, ExcessHet > 2)
hist(hets_highexcesshet$het_prop)
# 1 heterozygote would trigger this

hets_highexcesshet2 <- subset(hets, ExcessHet > 5)
head(hets_highexcesshet2[which(hets_highexcesshet2$het_prop < 0.2),])
hist(hets_highexcesshet2$het_prop)


### if I filtered on heterozygote proportion:
hets_filter1 <- subset(hets, het_prop < 0.2)
hist(hets_filter1$ExcessHet)
hist(hets_filter1[which(hets_filter1$ExcessHet < 20),"ExcessHet"])
hist(hets_filter1[which(hets_filter1$ExcessHet < 10),"ExcessHet"])


#### INBREEDING COEFFICIENT #####

hist(hets$InbreedingCoeff)

head(hets[which(hets$InbreedingCoeff > 0.5),])


#### MORE THAN 1 HETEROZYGOTE #####

head(hets)

hets_multiple <- subset(hets, HET > 1)
length(hets_multiple$CHROM) #2,642,278

hist(hets_multiple$ExcessHet)
hist(hets_multiple[which(hets_multiple$ExcessHet < 20),"ExcessHet"])
hist(hets_multiple[which(hets_multiple$ExcessHet < 10),"ExcessHet"])
head(hets_multiple[which(hets_multiple$ExcessHet < 10 &
                           hets_multiple$ExcessHet > 3 ),])

head(hets[which(hets$ExcessHet > 5),])

#### LESS THAN 20% MISSING #####

hets_missingfilter <- subset(hets, NCALLED > 230)
hist(hets_missingfilter$ExcessHet)
hist(hets_missingfilter[which(hets_missingfilter$ExcessHet < 20),"ExcessHet"])

# If 1 heterozygote out of 230, prop = 0.004
# ExcessHet would be 3.0197