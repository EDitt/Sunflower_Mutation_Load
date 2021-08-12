looking at the representation of the reference allele in the alignment

```bash
srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
module load R/4.0.0-foss-2019b
R
```


```R
dsnp <- read.table("Sunflower_SAM_Combined_Report.txt", sep = "\t", header=TRUE,
                     stringsAsFactors = FALSE)

dsnp$Alignment_list <- strsplit(dsnp$Alignment, "")

# using the "ReferenceAA" column from the combined report
dsnp$Refderived <- ifelse(apply(dsnp, 1, function(row) {
    row["ReferenceAA"] %in% unlist(row["Alignment_list"])
  }) == FALSE,
  "derived_state", "not_derived")

aggregate(dsnp$VariantID, by=list(dsnp$Refderived), length)
#        Group.1      x
# 1 derived_state      7
# 2   not_derived 645208

dsnp[which(dsnp$Refderived=="derived_state"),] # 2 have 'NA' for ReferenceAA, 5 have 'NA' for Alignment list


# frequency of reference in alignment across sites
dsnp$NumRefInAlignment <- apply(dsnp, 1, function(row) {
    length(which(unlist(row["Alignment_list"]) %in% row["ReferenceAA"]))
  }
)

dsnp$FreqRefInAlignment <- dsnp$NumRefInAlignment / dsnp$SeqCount

### are the ones where there is 1 represented in the alignment because of the reference?!?!
head(dsnp[which(dsnp$NumRefInAlignment==1),])


dsnp$pAdjusted <- p.adjust(dsnp$LogisticP_Masked, method = "BH", n = length(dsnp$LogisticP_Masked))

dsnp$Result <- ifelse(dsnp$pAdjusted < 0.05 & 
                      dsnp$SeqCount >= 10 & 
                      dsnp$MaskedConstraint < 1,
                      "Deleterious", "Tolerated")

# how many of the deleterious snps include only 1 in alignment?
aggregate(dsnp$VariantID, by=list(dsnp$Result), length)
#      Group.1      x
# 1 Deleterious 121877
# 2   Tolerated 523336

length(dsnp[which(dsnp$NumRefInAlignment==1),"VariantID"]) # 154,759

length(dsnp[which(dsnp$NumRefInAlignment==1 &
	dsnp$Result=="Deleterious"),"VariantID"]) # 11,796

# what proportion of 1's are deleterious?
aggregate(dsnp[which(dsnp$NumRefInAlignment==1), "VariantID"], by=list(dsnp[which(dsnp$NumRefInAlignment==1), "Result"]),
	length) # 8%
#      Group.1      x
# 1 Deleterious  11796
# 2   Tolerated 142963

aggregate(dsnp[which(dsnp$NumRefInAlignment>1), "VariantID"], by=list(dsnp[which(dsnp$NumRefInAlignment>1), "Result"]),
	length) # 22%
#      Group.1      x
# 1 Deleterious 110081
# 2   Tolerated 380368

##### SCRATCH

length(which(unlist(dsnp[1,"Alignment_list"]) %in% dsnp[1, "ReferenceAA"])) #76
length(grep(dsnp[1, "ReferenceAA"], unlist(dsnp[1,"Alignment_list"]), value=TRUE)) # 76 (another option)

length(unlist(dsnp[1, "Alignment_list"])) # = "SeqCount"

#####

```

Now how many do I get when I merge with VeP and check both the reference and alternate alleles?
(code from dSNP_table.R)
```R
vep <- read.table("/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm", sep = "\t", header=FALSE, stringsAsFactors = FALSE, na.strings = c("NA", "-"))

colnames(vep) <- c("VariantID", "Position", "Allele", "GeneID", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra")
mem_used() # 11.9 GB

vep_missense <- subset(vep, Consequence=="missense_variant")
rm("vep")

dsnp_data <- merge(dsnp, vep_missense, by=c("VariantID", "GeneID"))
dsnp_data$Amino_acidsSPLIT <- strsplit(dsnp_data$Amino_acids, "/") # split amino acids column
#dsnp_data$RefAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 1)) don't really need this because there's already a column for reference
dsnp_data$AltAA <- as.factor(sapply(dsnp_data$Amino_acidsSPLIT, "[", 2))
dsnp_data$Altderived <- ifelse(apply(dsnp_data, 1, function(row) {
    row["AltAA"] %in% unlist(row["Alignment_list"])
  }) == FALSE,
  "derived_state", "not_derived")
# this is new:
dsnp_data$Refderived <- ifelse(dsnp_data$NumRefInAlignment==1, "derived_state",
	ifelse(dsnp_data$NumRefInAlignment > 1, "not_derived",
		"NA"))
aggregate(dsnp_data$VariantID, by=list(dsnp_data$Refderived, dsnp_data$Altderived), length)
#        Group.1       Group.2      x
# 1 derived_state derived_state  87194
# 2            NA derived_state      5
# 3   not_derived derived_state 197740
# 4 derived_state   not_derived  67565
# 5   not_derived   not_derived 292709

# total alt derived (not including 5 NAs): 284,934

dsnp_data$Result <- ifelse(dsnp_data$pAdjusted < 0.05 & 
                               dsnp_data$SeqCount >= 10 & 
                               dsnp_data$MaskedConstraint < 1 & 
                               (dsnp_data$Refderived == "derived_state" | 
                                  dsnp_data$Altderived == "derived_state"),
                             "Deleterious", "Tolerated")

aggregate(dsnp_data$VariantID, by=list(dsnp_data$Result, dsnp_data$Refderived, dsnp_data$Altderived), length)
#      Group.1       Group.2       Group.3      x
#1   Tolerated derived_state derived_state  87194
#2   Tolerated            NA derived_state      5
#3 Deleterious   not_derived derived_state  76095
#4   Tolerated   not_derived derived_state 121645
#5 Deleterious derived_state   not_derived  11796
#6   Tolerated derived_state   not_derived  55769
#7   Tolerated   not_derived   not_derived 292709

aggregate(dsnp_data$VariantID, by=list(dsnp_data$Result, dsnp_data$Altderived), length)
# total alt derived + deleterious: # 76,095

### why did the number of Alt-derived deleterious mutations increase?
prevTable <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table", sep = "\t", header=T)
aggregate(prevTable$VariantID, by=list(prevTable$Refderived, prevTable$Altderived), length)
#        Group.1       Group.2      x
#1 derived_state derived_state      5
#2   not_derived derived_state 284934
#3   not_derived   not_derived 360274 (67,565 more)

aggregate(prevTable$VariantID, by=list(prevTable$Result, prevTable$Refderived, prevTable$Altderived), length)
#      Group.1       Group.2       Group.3      x
#1   Tolerated derived_state derived_state      5
#2 Deleterious   not_derived derived_state  54445
#3   Tolerated   not_derived derived_state 230489
#4   Tolerated   not_derived   not_derived 360274

# total alt derived (not including 5 NAs): 284,934
# total alt derived + deleterious: 54,445 ??????
```

Why did the number of alternate derived, deleterious alleles increase?

```R
#dsnp2 <- read.table("Sunflower_SAM_Combined_Report.txt", sep = "\t", header=TRUE,
#                     stringsAsFactors = FALSE)

Testdf <- TolvDel_sites("Sunflower_SAM_Combined_Report.txt", "/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm", 0.05, 10, 1, "LogisticP_Masked")
aggregate(Testdf$VariantID, by=list(Testdf$Result, Testdf$Refderived, Testdf$Altderived), length)
#      Group.1       Group.2       Group.3      x
# 1   Tolerated derived_state derived_state      5
# 2 Deleterious   not_derived derived_state  76095
# 3   Tolerated   not_derived derived_state 208839
# 4   Tolerated   not_derived   not_derived 360274

#### back to 76,095!
```

try using the full script
```bash
module load R/4.0.0-foss-2019b

Rscript /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/dSNP_table.R \
/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt \
/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm \
0.05 \
10 \
1 \
Masked \
/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Temp_test/dsnp_data_TEST.table

grep "Deleterious" dsnp_data_TEST.table | wc -l # 54,445 why???

## testing the script with a list statement
Rscript /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/dSNP_table.R \
/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt \
/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm \
0.05 \
10 \
1 \
Masked \
/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Temp_test/dsnp_data_TEST2.table

# need quotations?
Rscript "/home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/dSNP_table.R" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt" \
"/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm" \
"0.05" \
"10" \
"1" \
"Masked" \
"/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Temp_test/dsnp_data_TEST3.table"

```

why the difference?
```R
Testdf2 <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Temp_test/dsnp_data_TEST.table", 
	sep = "\t", header=TRUE, stringsAsFactors = FALSE)
aggregate(Testdf2$VariantID, by=list(Testdf2$Result, Testdf2$Refderived, Testdf2$Altderived), length)
#      Group.1       Group.2       Group.3      x
#1   Tolerated derived_state derived_state      5
#2 Deleterious   not_derived derived_state  54445
#3   Tolerated   not_derived derived_state 230489
#4   Tolerated   not_derived   not_derived 360274

Testdf2$Result2 <- ifelse(Testdf2$pAdjusted < 0.05 & 
                               Testdf2$SeqCount >= 10 & 
                               Testdf2$MaskedConstraint < 1 & 
                               (Testdf2$Refderived == "derived_state" | 
                                  Testdf2$Altderived == "derived_state"),
                             "Deleterious", "Tolerated")
aggregate(Testdf2$VariantID, by=list(Testdf2$Result2, Testdf2$Refderived, Testdf2$Altderived), length)
# this now has the correct number (76,095) of deleterious alleles???


Testdf2_diff <- Testdf2[which(Testdf2$Result!=Testdf2$Result2),]
min(Testdf2_diff$pAdjusted)
max(Testdf2_diff$pAdjusted) #  0.04965184
min(Testdf2_diff$SeqCount) # 36
max(Testdf2_diff$MaskedConstraint) # 9.222441e-05
which(Testdf2_diff$Altderived=="not_derived") #0

length(Testdf2$VariantID) # 645,213
length(unique(Testdf2$VariantID)) # 641,519

########
DelTest2 <- Testdf2[which(Testdf2$Result=="Deleterious"), "VariantID"]
DelTest1 <- Testdf[which(Testdf$Result=="Deleterious"), "VariantID"]
DelDiff <- setdiff(DelTest1, DelTest2)
length(DelDiff) # 21605

head(Testdf2[which(Testdf2$VariantID %in% DelDiff),])
head(Testdf[which(Testdf$VariantID %in% DelDiff),])
aggregate(Testdf[which(Testdf$VariantID %in% DelDiff),"VariantID"], 
	by=list(Testdf[which(Testdf$VariantID %in% DelDiff),"Result"]), length) #191 are tolerated, 21,618 are deleterious??
### duplicates

aggregate(Testdf2[which(Testdf2$VariantID %in% DelDiff),"VariantID"], 
	by=list(Testdf2[which(Testdf2$VariantID %in% DelDiff),"Result"]), length) # all tolerated (21,809)




# check p-values
# adjusted p-values are the same
Testdf2_diff <- Testdf2[which(Testdf2$VariantID %in% DelDiff),]
max(Testdf2_diff$pAdjusted) # 0.9999 ????
min(Testdf2_diff$pAdjusted) # 0.0261

Testdf1_diff <- Testdf[which(Testdf$VariantID %in% DelDiff & Testdf$Result=="Deleterious"),]
max(Testdf1_diff$pAdjusted) # 0.04965184
min(Testdf1_diff$pAdjusted) # 0.0261

#PvalDiff <- Testdf$pAdjusted - Testdf2$pAdjusted

both <- merge(Testdf[,c(1,7,13,31:34)], Testdf2[,c(1,7,13,31:34)], by="VariantID")

both$Pdiff <- both$pAdjusted.x - both$pAdjusted.y


head(both[which(both$Pdiff > 0.1 &
	both$VariantID %in% DelDiff),])

both$Pdiff <- apply(both, 1, function(row) {
    row["pAdjusted.y"] - row["pAjusted.x"]
    }
 )
both$Pdiff <- ifelse(both$pAdjusted.x != both$pAdjusted.y,
	both$pAdjusted.x - both$pAdjusted.y, 0)

Test3 <- read.table("/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Temp_test/dsnp_data_TEST2.table", 
	sep = "\t", header=TRUE, stringsAsFactors = FALSE)
```
everything looks ok-

'data.frame':   1 obs. of  6 variables:
 $ predict          : chr "/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt"
 $ vep_file         : chr "/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm"
 $ P.cutoff         : num 0.05
 $ mininumseq       : int 10
 $ maximumconstraint: int 1
 $ pvalcolumn       : chr "LogisticP_Masked"
