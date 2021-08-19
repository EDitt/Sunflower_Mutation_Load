# previous code that was in "Post_processing.md" but should probably go somewhere else



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

Use bcftools to count the number of alternate alleles per genotype
```bash
cd /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
module load BCFtools/1.10.2-GCC-8.3.0

bcftools stats -s - SAM_deleterious_polarized.vcf > DerivedDeleteriousStatsperSample.txt
grep "PSC" DerivedDeleteriousStatsperSample.txt > DerivedDeleteriousperSampleCounts.txt
```

# Polarize SNPs
I have a H. debilis ancestral FASTA sequence generated from ANGSD. Is it possible to use bcftools norm to fix the reference?

Ran `NormalizeVCF.sh` using the `--do-not-normalize` flag to just fix the reference (hopefully convert the reference to the ancestral allele?). I'm not sure what happens when it encounters an 'N' at a position. <- makes it multi-allelic with the 2 SNPs represented as 2 alternates (and literally 'N' is reference allele).
Submitted batch job 2692399

23,864,184 sites (out of 37,114,333) are multiallelic  the remaining (13,250,149) can be polarized by ancestral state
* the job failed due to issues with scaffold sequence so file is truncated, but will continue anyway (redo later)

### Filter for only biallelic positions

```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=12:00:00 --job-name=qlogin /bin/bash -l

bcftools view -m2 -M2 -v snps --threads 4 /scratch/eld72413/SAM_seq/dSNP_results/Sunflower_SAM_SNP_DebilisPolarized.vcf.gz --output-type v --output-file /scratch/eld72413/SAM_seq/dSNP_results/Sunflower_SAM_SNP_DebilisPolarized_BIALLELIC.vcf

grep -v "#" Sunflower_SAM_SNP_DebilisPolarized_BIALLELIC.vcf | wc -l # 13,250,148

# compress
sbatch --export=file=/scratch/eld72413/SAM_seq/dSNP_results/Sunflower_SAM_SNP_DebilisPolarized_BIALLELIC.vcf gzip_vcf.sh # Submitted batch job 2710927
```

### Subset this VCF file

Subset vcf file to make separate vcfs of both deleterious and tolerated
```bash
# deleterious
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Deleterious_positions.txt',vcf='/scratch/eld72413/SAM_seq/dSNP_results/Sunflower_SAM_SNP_DebilisPolarized_BIALLELIC.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/dSNP_results',name='SAM_deleterious_polarized' Subset_vcf.sh # Submitted batch job 2710964

grep -v "#" SAM_deleterious_polarized.vcf | wc -l # 39,010

# tolerated
sbatch --export=positions='/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Tolerated_positions_NoDups.txt',vcf='/scratch/eld72413/SAM_seq/dSNP_results/Sunflower_SAM_SNP_DebilisPolarized_BIALLELIC.vcf.gz',outputdir='/scratch/eld72413/SAM_seq/dSNP_results',name='SAM_tolerated_polarized' Subset_vcf.sh # Submitted batch job 2710965

grep -v "#" SAM_tolerated_polarized.vcf | wc -l # 398,104
```

### Frequency distributions
```bash
module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
vcftools --gzvcf /scratch/eld72413/SAM_seq/dSNP_results/Sunflower_SAM_SNP_DebilisPolarized_BIALLELIC.vcf.gz --freq --positions /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Deleterious_positions.txt --out SAM_deleterious_polarized

vcftools --gzvcf /scratch/eld72413/SAM_seq/dSNP_results/Sunflower_SAM_SNP_DebilisPolarized_BIALLELIC.vcf.gz --freq --positions /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Tolerated_positions_NoDups.txt --out SAM_tolerated_polarized

module load R/4.0.0-foss-2019b
R
```

```R
library(tidyr)
hist_breaks <- seq(0, 1, by=0.05)

del <- read.table("SAM_deleterious_polarized.frq", sep = "\t", header=TRUE, row.names=NULL)
colnames(del) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "Ancestral_Freq", "Derived_Freq")
del <- del %>% separate(Derived_Freq, c("Derived_Base", "Derived_Freq"), sep=":")
del$Derived_Freq <- as.numeric(del$Derived_Freq)

SAM_deleterious_POLARIZED_freqbins <- hist(del$Derived_Freq, plot=FALSE, breaks=hist_breaks)
save(SAM_deleterious_POLARIZED_freqbins, file="DeleteriousBin_POLARIZED.RData")

tol <- read.table("SAM_tolerated_polarized.frq", sep = "\t", header=TRUE, row.names=NULL)
colnames(tol) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "Ancestral_Freq", "Derived_Freq")
tol <- tol %>% separate(Derived_Freq, c("Derived_Base", "Derived_Freq"), sep=":")
tol$Derived_Freq <- as.numeric(tol$Derived_Freq)

SAM_tolerated_POLARIZED_freqbins <- hist(tol$Derived_Freq, plot=FALSE, breaks=hist_breaks)
save(SAM_tolerated_POLARIZED_freqbins, file="ToleratedBin_POLARIZED.RData")

# on local computer
```

Use bcftools to count the number of alternate alleles per genotype

```bash
module load BCFtools/1.10.2-GCC-8.3.0

bcftools stats -s - SAM_deleterious_polarized.vcf > DerivedDeleteriousStatsperSample.txt
grep "PSC" DerivedDeleteriousStatsperSample.txt > DerivedDeleteriousperSampleCounts.txt
```
### Polarize SNPs using ANGSD-

I tried installing angsd-
```bash
module load  HTSlib/1.10.2-GCC-8.3.0
wget http://popgen.dk/software/download/angsd/angsd0.934.tar.gz
tar xf angsd0.934.tar.gz
cd angsd
make
```
there is a now this exectuable: `/home/eld72413/apps/angsd/angsd`
I can use this to run Chaochih's code to make an ancestral vcf file: https://github.com/ChaochihL/Barley_Outgroups/blob/master/morex_v1/angsd_anc_inf.job



########################################################
### Previous numbers before fixing Predict output parsing


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