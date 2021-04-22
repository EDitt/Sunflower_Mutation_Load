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
