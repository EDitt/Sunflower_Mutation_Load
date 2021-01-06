setwd("/Users/eld72413/Google Drive/Active Projects/DelMutation/Variant_graphs")

hets <- read.csv("Variants_HetInfo.table", sep = "\t", header = TRUE)
length(hets$CHROM) #4,903,270
hets$het_prop <- hets$HET / hets$NCALLED

hist(hets$het_prop)

length(which(hets$het_prop > 0.2)) # 166,638 (3.4%)

hist(hets$ExcessHet)
hist(hets[which(hets$ExcessHet < 20),"ExcessHet"])
hist(hets[which(hets$ExcessHet < 10),"ExcessHet"])

HighHetProp <- which(hets$het_prop > 0.2)
HighExcessHet <- which(hets$ExcessHet > 5)
length(HighExcessHet) #395,838 (8%)
length(which(HighHetProp %in% HighExcessHet)) #112,639
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