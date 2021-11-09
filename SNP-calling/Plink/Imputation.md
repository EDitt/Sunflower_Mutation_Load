# Imputation for GWAS

Tutorial: https://github.com/adrianodemarino/Imputation_beagle_tutorial

### Software Needed

```bash
# software needed:
module load BCFtools/1.13-GCC-8.3.0
module load R/4.1.0-foss-2019b
module load Eagle/2.4.1-linux_x86_64
module load BEAGLE/5.2-Java-1.8.0_144
```

Check whether R packages are availabile:
```R
library(data.table) # data.table 1.14.2 using 8 threads (see ?getDTthreads).  Latest news: r-datatable.com
library(sm) # Package 'sm', version 2.2-5.7: type help(sm) for summary information
```

### Files Needed

Reference Genome & Genetic Map

#### Reference Genome Files
```bash
GenomeFASTA=/scratch/eld72413/SunflowerGenome/Ha412HOv2.0-20181130.fasta
GenomeINDEX=Ha412HOv2.0-20181130.fasta/Ha412HOv2.0-20181130.fasta.fai

```

#### Genetic Map Files

Genetic Map File Format?
According to: https://mathgen.stats.ox.ac.uk/impute/input_file_options.html
should have 3 columns: physical position (in bp), recombination rate between current position and next position in map (in cM/Mb), and genetic map position (cM). Header line with an unbroken character string for each column (e.g., "position COMBINED_rate(cM/Mb) Genetic_Map(cM)")

###### Genetic Map File for Phasing
Downloaded the example file: genetic_map_hg38_withX.txt.gz
```bash
gunzip genetic_map_hg38_withX.txt.gz 
head genetic_map_hg38_withX.txt # chr, position, COMBINED_rate(cM/Mb), Genetic_Map(cM)
# chromosome is labelled "1-23"

# modified this text from tutorial (since I unzipped the file):
for CHR in {1..23}; do
    cat genetic_map_hg38_withX.txt | \
    grep ^${CHR} | \
    sed '1ichr position COMBINED_rate(cM/Mb) Genetic_Map(cM)' \
    > eagle_chr${CHR}_b38.map
done
### this just split the file across chromosomes
```

###### Genetic Map File for Imputation <- I think this is the raw variant file that imputation will be performed on?
```bash
tar -x plink.GRCh38.map.zip
```

#### Reference Panel Files
vcf.gz and vcf.gz.tbi files
```bash
# example file
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{1..22},X}_GRCh38.genotypes.20170504.vcf.gz{,.tbi} # not found
```

Looked at plink.imiss file to see the distribution of missing data across genotypes
Subset for the genotypes with the least missing data
```R
setwd("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Plink")
plink_geno <- read.table("plink.imiss", header = TRUE)
stats <- boxplot(plink_geno$F_MISS, plot=FALSE) # lower hinge of boxplot: 0.033150
length(plink_geno[which(plink_geno$F_MISS <= 0.033150),"IID"]) # 72

GoldStandard <- plink_geno[which(plink_geno$F_MISS <= 0.033150), c(1,2)]

write.table(GoldStandard, file = "LowMissingGeno.txt", 
            quote=FALSE, 
            row.names = FALSE, col.names = FALSE)

### to filter with bcftools
GoldStandard2 <- as.character(plink_geno[which(plink_geno$F_MISS <= 0.033150), 2])
# need to change the name of 2 genotypes:
GoldStandard2[which(GoldStandard2=="33")] <- "SF_33"
GoldStandard2[which(GoldStandard2=="PPN136")] <- "NMS373_PPN136"

write.table(GoldStandard2, file = "LowMissingGeno_bcftools.txt", 
            quote=FALSE, 
            row.names = FALSE, col.names = FALSE)
```

Filter VCF file for these genotypes
```bash
module load PLINK/1.9b_5-x86_64
plink --file /scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2_FILTERED \
--keep /scratch/eld72413/SAM_seq/Plink/LowMissingGeno.txt \
--missing --allow-extra-chr

awk '{if ($3==0) {print $0}}' plink.lmiss | wc -l #3,948,880

# Make new filtered file
plink --file /scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2_FILTERED \
--keep /scratch/eld72413/SAM_seq/Plink/LowMissingGeno.txt \
--make-bed --allow-extra-chr

```

261 set of genotypes in salt data
```bash
sed -i 's/SAM/PPN/g' Genotypes_261set.txt

plink --file /scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2_FILTERED \
--keep /scratch/eld72413/SAM_seq/Plink/Genotypes_261set.txt \
--missing --allow-extra-chr

awk '{if ($3==0) {print $0}}' plink.lmiss | wc -l #248
```

### Pre-processing for Reference Pannel
- rename numerical chromosome names with 'chr' tag
- remove rare variants (AC threshold with bcftools)
- split multiallelic sites to biallelic records (bcftools norm)
- Keep only SNPs and INDELs (exclude structural variants)
- Align variants to ref genome (bcftools norm)
- Remove duplicate variants
- Remove multiallelic records (after bcftools norm)
- Remove sites containing missing data

I have already done all of these (not sure about the chromosomal naming) except for removing missing data
I will first, subset the VCF to 72 genotypes with the least missing data (< 3.315% (lower hinge of boxplot) of sites missing)
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l

module load BCFtools/1.10.2-GCC-8.3.0

bcftools view -S /scratch/eld72413/SAM_seq/Plink/LowMissingGeno_bcftools.txt \
/scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2_missing_filtered.vcf \
-Ou | \
bcftools view -g ^miss -Oz -o /scratch/eld72413/SAM_seq/Plink/Sunflower_SAM_HA412v2_REF_PANEL.vcf.gz

bcftools stats Sunflower_SAM_HA412v2_REF_PANEL.vcf.gz # 3,948,880 SNPs

####
#bcftools filter -i 'ID=@/scratch/eld72413/SAM_seq/Plink/LowMissingGeno.txt' \




```

### Reference panel allele frequencies:
Generate a tab-delimited file of the reference panel allele frequencies, one variant per line, with columns CHR, SNP (in format CHR_POS_REF_ALT), REF, ALT, AF (including the header line)

### Create binary reference panel files
the phased reference panel files per chromosome are required in bref format

