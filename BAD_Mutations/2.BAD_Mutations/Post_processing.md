
# Filter VeP Report

```bash
module load VEP/101.0-foss-2019b-Perl-5.30.0
INPUT=/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm
filter_vep -i ${INPUT} -o fullsam_missense.txt -filter "Consequence is missense_variant"
filter_vep -i ${INPUT} -o fullsam_synon.txt -filter "Consequence is synonymous_variant"

# remove header lines to use with R
#### actually R automatically ignores lines that start with "#"
awk 'NR > 29 {print}' fullsam_missense.txt > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/fullsam_missense_noHEADER.txt
awk 'NR > 29 {print}' fullsam_synon.txt > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/fullsam_synon_noHEADER.txt
```

# Significant Variants from Compiled Predict report

I tested 50,838 codons for Sunflower, so significance threshold= 0.05/50838 <- way too stringent for this dataset

Will also filter out alignments with fewer than 10 species

I ended up writing an R script to do what I was trying to do

```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results

srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
# job killed at mem=22g
module load R/4.0.0-foss-2019b
module load R_ML/3.3.3 # for MSI

Rscript /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/2.BAD_Mutations/dSNP_table.R \
/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Sunflower_SAM_Combined_Report.txt \
/scratch/eld72413/SAM_seq/VeP/SAM_SNP_Final_BiallelicNorm \
0.05 \
10 \
1 \
Masked \
/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/dsnp_data.table

```

How many are deleterious vs. tolerated?
Get positions to subset VCF
```bash
grep "Tolerated" dsnp_data.table | wc -l # 590,768
grep "Deleterious" dsnp_data.table | wc -l # 54,445

## deleterious positions
grep "Deleterious" dsnp_data.table | awk '{print $17}' | awk '{$1=$1}1' FS=':' OFS='\t' > Deleterious_positions.txt
## tolerated positions
grep "Tolerated" dsnp_data.table | awk '{print $17}' | awk '{$1=$1}1' FS=':' OFS='\t' > Tolerated_positions.txt

cd /scratch/eld72413/SAM_seq/VeP
grep -v "#" fullsam_synon.txt | awk '{print $2}' | awk '{$1=$1}1' FS=':' OFS='\t' > /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positions.txt
```

# Remove duplicate positions

Some variant ID's are the same among the list of positions if they code for different genes on the forward and reverse strands. I will use the most severe consequence of the positions

Any duplicates within a category will automatically be removed when I subset the VCF by position

```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load R/4.0.0-foss-2019b
```
Use R
```R
synon <- read.table("Synonymous_positions.txt", sep = "\t", header=FALSE,
                         stringsAsFactors = FALSE)
synon$pos <- paste0(synon$V1,"_",synon$V2)

tolerated <- read.table("Tolerated_positions.txt", sep = "\t", header=FALSE,
                         stringsAsFactors = FALSE)
tolerated$pos <- paste0(tolerated$V1,"_",tolerated$V2)

dsnp <- read.table("Deleterious_positions.txt", sep = "\t", header=FALSE,
                         stringsAsFactors = FALSE)
dsnp$pos <- paste0(dsnp$V1,"_",dsnp$V2)

# which dsnp positions are also considered tolerated
length(which(tolerated$pos %in% dsnp$pos)) #525
tolerated_nodups <- tolerated[-c(which(tolerated$pos %in% dsnp$pos)),]
length(tolerated$pos) # 590,768
length(tolerated_nodups$pos) # 590,243
write.table(tolerated_nodups[,c(1:2)], "Tolerated_positions_NoDups.txt", sep = "\t",
	col.names=FALSE, row.names=FALSE, quote=FALSE)

# which synonymous positions are considered tolerated OR deleterious?
length(which(synon$pos %in% tolerated$pos)) # 4593
length(which(synon$pos %in% dsnp$pos)) #396
length(which(synon$pos %in% tolerated$pos |
		synon$pos %in% dsnp$pos)) # 4989

synon_nodups <- synon[-c(which(synon$pos %in% tolerated$pos |
							synon$pos %in% dsnp$pos)),]
length(synon$pos) # 835,622
length(synon_nodups$pos) #830,633

write.table(synon_nodups[,c(1:2)], "Synonymous_positions_NoDups.txt", sep = "\t",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
```
*** I will need to remove duplicate positions for syonymous variants (not also annotated as missense)
Write a script to parse VeP output to remove duplicate positions (go by most severe annotation)

Check
```bash
wc -l Synonymous_positions_NoDups.txt # 830,633
awk '{print $0}' Synonymous_positions_NoDups.txt | sort -u | wc -l # 827,101 <- number of positions to expect in VCF
awk '{print $0}' Synonymous_positions_NoDups.txt | sort | uniq -cd | wc -l # 3532 are duplicates

wc -l Tolerated_positions_NoDups.txt # 590,243
awk '{print $0}' Tolerated_positions_NoDups.txt | sort -u | wc -l # 587,108 <- number of positions to expect in VCF
awk '{print $0}' Tolerated_positions_NoDups.txt  | sort | uniq -cd | wc -l # 3135 are duplicates

wc -l Deleterious_positions.txt # 54,445
awk '{print $0}' Deleterious_positions.txt | sort -u | wc -l # 54,411 <- number of positions to expect in VCF
awk '{print $0}' Deleterious_positions.txt  | sort | uniq -cd | wc -l #34 are duplicates

wc -l Sunflower_SAM_Combined_Report.txt # 645,216 (includes header line)
```

Total Number of unique positions= 
Synonymous: 827,101

Tolerated: 587,108
Deleterious: 54,411
		Total (unique) Missense predicted: 641,519

Number in Predict output: 645,215 (difference= 3,696)
	- 3135 duplicate positions (tolerated)
	- 34 duplicate positions (deleterious)
	- 525 positions annotated as deleterious & tolerated
	Total: 3,694 (there are 2 NA's in predict output)

# Subset VCF file

Subset vcf file to make separate vcfs of both deleterious and tolerated
```bash
# deleterious
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Deleterious_positions.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_deleterious' Subset_vcf.sh # Submitted batch job 2245897

grep -v "#" SAM_deleterious.vcf | wc -l # 54,411

# tolerated (redo to remove duplicate positions)
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Tolerated_positions_NoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_tolerated' Subset_vcf.sh # Submitted batch job 2378842

grep -v "#" SAM_tolerated.vcf | wc -l # 587,108  (before removing positions in dSNPs: 587,633)
```

### Synonymous positions
Subset VCF to get synonymous positions

```bash
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Synonymous_positions_NoDups.txt',vcf='/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/VarFilter_All/Sunflower_SAM_SNP_Calling_BIALLELIC_norm.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results',name='SAM_synonymous' Subset_vcf.sh # 2378848


grep -v "#" SAM_missense.vcf | wc -l 
wc -l Synonymous_positions.txt
```



