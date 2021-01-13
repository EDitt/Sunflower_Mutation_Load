## Variant Analysis

I ran VeP on the full set (including multi-allelic sites). 
- Used `VeP.sh` script

Calculate site frequency spectrum with vcftools
- Used `Freq.sh`
- ~~Took output from that and graphed in R~~ file too big
```bash
cd /scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221

tmux new -s Rcode
srun --pty  -p inter_p  --mem=4G --nodes=1 --ntasks-per-node=1 --time=6:00:00 --job-name=qlogin /bin/bash -l #1027971

module load R/4.0.0-foss-2019b
module load Python/3.8.2-GCCcore-8.3.0

# use scripts from sequence handling
seqhand="/home/eld72413/MorrellSeqHandling/sequence_handling"
vcf="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Sunflower_SAM_SNP_Calling_Final_Filtered.vcf"
out="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Stats"
name="Sunflower_SAM_SNP_Calling_Final"

python3 "${seqhand}/HelperScripts/VCF_MAF.py" "${vcf}" > "${out}/${name}_MAF.txt"

# this did not complete
srun: Force Terminated job 1027971
srun: error: c1-7: task 0: Out Of Memory
```

Try bcftools stats plotting functions
```bash
module load BCFtools/1.10.2-GCC-8.3.0
module load matplotlib/3.1.1-intel-2019b-Python-3.7.4 #need for this function
plot-vcfstats -p ${out} Sunflower_SAM_SNP_FinalStats.txt

# error:
Neither pdflatex or tectonic were found in your PATH, impossible to create a PDF at /apps/eb/BCFtools/1.10.2-GCC-8.3.0/bin/plot-vcfstats line 111.
        main::error('Neither pdflatex or tectonic were found in your PATH, impossi...') called at /apps/eb/BCFtools/1.10.2-GCC-8.3.0/bin/plot-vcfstats line 1820
        main::create_pdf('HASH(0x25b5ad8)') called at /apps/eb/BCFtools/1.10.2-GCC-8.3.0/bin/plot-vcfstats line 72
# can try editing the script it put in the output
```

Summarize substitution type information from bcftools stats
```R
library(ggplot2)

Substitutionsdf <- data.frame(
  Substitutions = c("[A/T]", "[A/C]=[T/G]", "[A/G]=[T/C]", "[C/G]"),
  Number = c(6418253, 10567791, 35033872, 3423934)
  )
Substitutionsdf$Substitutions <- factor(Substitutionsdf$Substitutions, 
                                        levels=c("[A/T]", "[A/C]=[T/G]", "[A/G]=[T/C]", "[C/G]"))

ggplot(Substitutionsdf, aes(x=Substitutions, y=Number)) + 
  geom_bar(stat = "identity")
Substitutionsdf$Proportion <- Substitutionsdf$Number / sum(Substitutionsdf$Number)
ggplot(Substitutionsdf, aes(x=Substitutions, y=Proportion)) + 
  geom_bar(stat = "identity") +
  theme_bw()

```


Compress the file (including multi-allelic sites)
```bash
# tmux window: gzip (?)
module load BCFtools/1.10.2-GCC-8.3.0
bcftools view Sunflower_SAM_SNP_Calling_Final_Filtered.vcf -Oz -o Sunflower_SAM_SNP_Calling_Final_Filtered.vcf.gz
# *compression completed before interactive job ran out of walltime
```

Will look at concordance with truth set
```bash
tmux new -s concordance
module load BCFtools/1.10.2-GCC-8.3.0

# first need to compress + index truth set
Truth_Set="/scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf"

bcftools view $Truth_Set -Oz -o /scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf.gz
bcftools index /scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf.gz

Truth_zip="/scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf.gz"

srun --pty  -p inter_p  --mem=2G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l

VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Sunflower_SAM_SNP_Calling_Final_Filtered.vcf.gz"
Truth_zip="/scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf.gz"
OutputDir="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Stats"

# need to index gzipped vcf file
bcftools index $VCF

bcftools stats $VCF $Truth_zip > ${OutputDir}/ConcordanceStats_SNParray6524.txt
```

Number of SNPs only in my set: 54,140,752
Number of SNPs only in truth set: 1,496
Number of shared SNPs: 5,028

Recovered 77% of SNPs in array
I should figure out what the annotation values are for the 1,496 SNPs that aren't recovered
- Ran script `gatk_SelectConcordant.sh`

## Biallelic sites only

Filter to biallelic sites
```bash
tmux new -s biallelic
srun --pty  -p inter_p  --mem=8G --nodes=1 --ntasks-per-node=4 --time=24:00:00 --job-name=qlogin /bin/bash -l #1027883

module load BCFtools/1.10.2-GCC-8.3.0
OUTPUT_DIR="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Biallelic"
VCF="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter6_011221/Sunflower_SAM_SNP_Calling_Final_Filtered.vcf"

bcftools view -m2 -M2 -v snps --threads 4 ${VCF} --output-type v --output-file ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_BIALLELIC.vcf

# this completed
# how many sites?
bcftools stats ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_BIALLELIC.vcf > BiallelicSTATS.txt
# 51,014,412, ts/tv ratio is 1.82

# compress biallelic site vcf file *** currently in progress
bcftools view ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_BIALLELIC.vcf -Oz -o ${OUTPUT_DIR}/Sunflower_SAM_SNP_Calling_BIALLELIC.vcf.gz

# look at concordance with truth set
Truth_zip="/scratch/eld72413/SNParray/FinalFiles/MapUniqueSNP_idt90_rename_rmContigs_sorted.vcf.gz"

```


