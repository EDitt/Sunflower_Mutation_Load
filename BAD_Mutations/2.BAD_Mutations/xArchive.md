# Code that didn't work/needs troubleshooting/no longer needed- all related to downstream analyses

### Add SNP name

Add SNP name to annotated variants using the Substitution_ID.py script in the 'Supporting' directory of BAD_Mutations
```bash
module load python3/3.6.3_anaconda5.0.1
cd /panfs/roc/groups/9/morrellp/shared/Software/BAD_Mutations/Supporting
compiled=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Compile_58perc/Sunflower_SAM_Combined_Report.txt
sub_file=/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/SAM_SNP_BadMut_Summary_edit
python Substitution_ID.py $compiled $sub_file
```
Traceback (most recent call last):
  File "Substitution_ID.py", line 87, in <module>
    tmp[0] = subs[key]
KeyError: ('Ha412HOChr00c00005g0858651', '63')

```bash
grep Ha412HOChr00c00005g0858651 $sub_file
```
Ha412HOChr00c00005g0858651	63	R	Ha412HOChr00c00005_182708_A/G

### Put VeP results in table

Cloned Deleterious_Mutations repository: https://github.com/MorrellLAB/Deleterious_Mutations
```bash
module load Biopython/1.75-foss-2019b-Python-2.7.16
Del_Mut_Path="/home/eld72413/DelMut/Deleterious_Mutations/Analysis_Scripts"
vcf=/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz
vep=/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm
python ${Del_Mut_Path}/VEP_To_Table.py ${vcf} ${vep} > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/VEP_Table

```
^ this did not work, possibly only meant for soybean samples?
Traceback (most recent call last):
  File "/home/eld72413/DelMut/Deleterious_Mutations/Analysis_Scripts/VEP_To_Table.py", line 115, in <module>
    main()
  File "/home/eld72413/DelMut/Deleterious_Mutations/Analysis_Scripts/VEP_To_Table.py", line 90, in main
    vdata = parse_vcf(sys.argv[1])
  File "/home/eld72413/DelMut/Deleterious_Mutations/Analysis_Scripts/VEP_To_Table.py", line 19, in parse_vcf
    snpid = tmp[2]
IndexError: list index out of range

### R code I used before developing script

```R
dsnp <- read.table("Sunflower_SAM_Combined_Report_preliminary2.txt", sep = "\t", header=TRUE)

missense <- read.table("fullsam_missense_noHEADER.txt", sep = "\t", header=FALSE)
colnames(missense) <- c("VariantID", "Position", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
# columns I need:
# 1) variant ID that matches with the BAD_Mutations report
# 2) snp position
# 3) variant allele
# 4) gene
# 11) amino acids

# merge based on position
dsnp_data <- merge(dsnp, missense, by="VariantID")

# split the Amino_acids column:
dsnp_data$Amino_acidsSPLIT <- strsplit(dsnp_data$Amino_acids, "/")
dsnp_data$RefAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 1))
dsnp_data$AltAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 2))

# I think(?) the "ReferenceAA" in the BAD_Mutations output is the reference of unrelated Angiosperm genomes
length(which(dsnp_data$RefAA != dsnp_data$ReferenceAA | dsnp_data$AltAA != dsnp_data$ReferenceAA)) # 425,939
length(which(dsnp_data$RefAA != dsnp_data$ReferenceAA | dsnp_data$AltAA != dsnp_data$ReferenceAA &
	dsnp_data$MaskedConstraint < 1 )) # 350,284

# Tom's criteria
lrt_sig <- 0.05/50838
length(which(dsnp_data$LogisticP_Masked < lrt_sig)) # 0
min(na.omit(dsnp_data$LogisticP_Masked)) # 3.345465e-05
min(na.omit(dsnp_data$LogisticP_Unmasked)) # 3.310914e-07

# too stringent, try FDR
dsnp_data$pAdjusted <- p.adjust(dsnp_data$LogisticP_Masked, method = "BH", n = length(dsnp_data$LogisticP_Masked))
length(which(dsnp_data$pAdjusted < 0.05)) # 83642

minseq <- 10
length(which(dsnp_data$SeqCount >= minseq)) # 385827 (out of 420,768)

length(which(dsnp_data$SeqCount >= minseq & dsnp_data$pAdjusted < 0.05)) # 83642

max_constraint <- 1

lrt <- dsnp_data[(dsnp_data$pAdjusted < 0.05 & dsnp_data$SeqCount >= minseq & 
	dsnp_data$MaskedConstraint < max_constraint & (dsnp_data$RefAA != dsnp_data$ReferenceAA | dsnp_data$AltAA != dsnp_data$ReferenceAA)), ] 
length(lrt$VariantID) # 81565

#lrt <- dsnp_data[(dsnp_data$pAdjusted < 0.05 & dsnp_data$SeqCount >= minseq & 
#	dsnp_data$MaskedConstraint < max_constraint), ] # 81565

#lrt <- dsnp_data[(dsnp_data$pAdjusted < 0.05 & dsnp_data$SeqCount >= minseq), ] # 83642

write.table(lrt[,-30], "dsnp_data_PRELIM.table", sep = "\t", quote=FALSE, row.names=FALSE)
```

# Using function in R 
```R
source("/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/dSNP_table.R")

Mydf <- TolvDel_sites("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report_preliminary2.txt",
"/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm",
0.05,
10,
1,
"LogisticP_Masked")

write.table(Mydf[-30], "dsnp_data_PRELIM_NEW.table", sep = "\t", quote=FALSE, row.names=FALSE)
```

# Figuring out issue with R script that was causing duplicates after merging with VeP output (and why there were duplicates in the VeP output)

```bash

wc -l SAM_SNP_BadMut_Summary # 704,075
awk '{print $4}' SAM_SNP_BadMut_Summary | sort -u | wc -l #699,805 # unique variants
awk '{print $1}' SAM_SNP_BadMut_Summary | sort -u | wc -l #50,838 # unique regions

wc -l dsnp_data.table # 1,120,517  ## why the difference?
awk 'NR > 1 {print $1}' dsnp_data.table | sort -u | wc -l # 641,519

wc -l Sunflower_SAM_Combined_Report.txt # 645,216 (including header line)

grep -v "#" SAM_SNP_Final_BiallelicNorm | wc -l # 43,019,114
grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{print $1}' | sort -u | wc -l # 37,120,112

grep "missense_variant" SAM_SNP_Final_BiallelicNorm | wc -l # 712,631
grep "missense_variant" SAM_SNP_Final_BiallelicNorm | awk '{print $1}' | sort -u | wc -l # 708,311 more than expected..
grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="missense_variant") {print $1}}' | wc -l # 704,075
grep -v "#" SAM_SNP_Final_BiallelicNorm | awk '{if ($7=="missense_variant") {print $1}}' | sort -u | wc -l # 699,805

grep -v "#" fullsam_missense.txt | wc -l # 704,075
grep -v "#" fullsam_missense.txt | awk '{print $1}' | sort -u | wc -l # 699,805
grep -v "#" fullsam_missense.txt | awk '{print $1}' | sort | uniq -cd | head
grep -v "#" fullsam_missense.txt | awk '{print $1}' | sort | uniq -cd | wc -l # 4270 are duplicated

```
Summary: there should be 699,805 unique missense variant positions in the VeP output (not counting some that have multiple consequences including missense).

Need to figure out where there are so many more positions in the dsnp_data.table
```R
missense <- read.table("/scratch/eld72413/SAM_seq/VeP/fullsam_missense.txt", sep = "\t", header=FALSE)
synon <- read.table("/scratch/eld72413/SAM_seq/VeP/fullsam_synon.txt", sep = "\t", header=FALSE)
colnames(missense) <- c("VariantID", "Position", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
colnames(synon) <- c("VariantID", "Position", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")

length(synon$VariantID) # 835,622
length(missense$VariantID) # 704,075

dsnp_table <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table", sep = "\t", header=TRUE)
length(dsnp_table$VariantID) # 1,120,516

length(which(dsnp_table$VariantID %in% missense$VariantID)) # 1,120,516
length(which(dsnp_table$VariantID %in% synon$VariantID)) # 12,962
length(which(missense$VariantID %in% synon$VariantID)) # 5346

dsnp <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt", sep = "\t", header=TRUE,
                     stringsAsFactors = FALSE)
length(dsnp$VariantID) # 645,215

length(which(dsnp$VariantID %in% missense$VariantID)) # 645213
dsnp[which(!dsnp$VariantID %in% missense$VariantID),] # 2 variant ID's are "NA"
length(which(dsnp$VariantID %in% synon$VariantID))

length(which(is.na(dsnp$VariantID))) # 2
length(which(is.na(dsnp$LogisticP_Masked)))
dsnp[which(is.na(dsnp$LogisticP_Masked)),]
dsnp[which(is.na(dsnp$VariantID)),]

length(which(is.na(dsnp_table$VariantID))) #0

# full vep output
vep <- read.table("/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm", sep = "\t", header=FALSE,
                         stringsAsFactors = FALSE, na.strings = c("NA", "-"))
colnames(vep) <- c("VariantID", "Position", "Allele", "GeneID", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
length(vep$VariantID) # 43,019,114
length(which(is.na(vep$VariantID))) #0
length(unique(vep$VariantID)) # 37,120,112

#vep$fConsequence <- as.factor(vep$Consequence)
#levels(vep$fConsequence) # there is one category that has 2 consequences (including missense)

vep_missense <- subset(vep, Consequence=="missense_variant")
length(vep_missense$VariantID) # 704,075
length(unique(vep_missense$VariantID) # 699,805

# try merging again
dsnp_data1 <- merge(dsnp, vep, by="VariantID")
length(dsnp_data1$VariantID) # 1,120,516

#dsnp_data2 <- merge(dsnp, vep[unique(vep$VariantID),], by="VariantID") 
#length(dsnp_data2$VariantID)

# merge using this code?
dsnp_data2 <- merge(dsnp, vep_missense, by=c("VariantID", "GeneID"))
length(dsnp_data2$VariantID)

```
In my compiled prediction report, there are 2 variant ID's that are listed as "NA".
The total number of compiled variants in this report is: 645,215 (including the 2 NA's). This is out of 704,075 missense variants in the VeP output. The only 2 that are missing are the "NA's".

For some reason, my R code gives me a table with 1,120,516 variants. These are all listed as being in the missense vep output (which only has 704,075 rows). Therefore, there must be duplicates.

To add to R script
```R
vep_missense <- subset(vep, Consequence=="missense_variant")
rm("vep")

dsnp_data <- merge(dsnp, vep_missense, by=c("VariantID", "GeneID"))
length(dsnp_data$VariantID) # 645,213 (this probably excludes the NA's)

#vep_missenseU <- vep_missense[unique(vep_missense$VariantID),] # this took a loooooong time. did not work (everything is NA)
#length(vep_missenseU$VariantID) # 699,805

library(pryr)
mem_used() # with vep, vep_missense, and dsnp: 11.5 GB
```

Duplicates are caused by different genes on forward and reverse strand. Example:

Ha412HOChr01_316258_A/T Ha412HOChr01:316258 T Ha412HOChr01g0000011  mRNA:Ha412HOChr01g0000011 Transcript  missense_variant  850 850 284 F/I Ttt/Att - IMPACT=MODERATE;STRAND=-1;SOURCE=Ha412HOv2.0.gff3.gz

Ha412HOChr01_316258_A/T Ha412HOChr01:316258 T Ha412HOChr01g0000021  mRNA:Ha412HOChr01g0000021 Transcript  missense_variant  239 6 2 K/N aaA/aaT - IMPACT=MODERATE;STRAND=1;SOURCE=Ha412HOv2.0.gff3.gz

Solution: need to also sort by Gene ID


# Initially subset VCF by missense positions (before I realized the duplicate position issue)

### Missense positions

First subset vcf to find *all* missense positions
```bash

awk 'NR > 1 {print $2}' fullsam_missense_noHEADER.txt | awk '{$1=$1}1' FS=':' OFS='\t' > Missense_positions.txt

sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Missense_positions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_missense' Subset_vcf.sh # 2164550


grep -v "#" SAM_missense.vcf | wc -l # 699,805
wc -l Missense_positions.txt # 704,075 # not sure why they aren't equal
```
