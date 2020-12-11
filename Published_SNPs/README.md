### Received raw version of SNP set from Greg Owens
Brian put in directory `/scratch/bp26123/UBC_snps`

I want to look at basic stats, SFS, and run VeP

```bash
SNPs="/scratch/bp26123/UBC_snps/Annuus.tranche90.snp.remappedHa412.vcf.gz"
```

Made a script to run bcftools (need to also test out SLURM as the queuing system was recently changed)
Number of samples: 1293
Number of records: 20,948,148
ts/tv: 1.32
number of singletons: 2,003,528

### Subset VCF by group (wild vs. SAM lines vs landrace)

First, subset for SAM lines

```bash
module load BCFtools/1.10.2-GCC-8.3.0

bcftools query -l $SNPs > SampleNames.txt #full set of samples in this SNP set (includes wild Annuus, SAM lines, and landraces)

# Used my VCF file to get sample names of SAM lines
SAM_SNPs_new="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/New2/Filter1_102120/Sunflower_SAM_SNP_Calling_snps.filtered.vcf"
bcftools query -l $SAM_SNPs_new > SAM_SampleNames288.txt

# need to change SAM to PPN to match 
sed -i 's/PPN/SAM/g' SAM_SampleNames288.txt
# changed RHA415-4_SAM251, Hopi_SAM285, NMS373_SAM136 manually
```

Subset VCF
```bash
module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

OUTPUTDIR="/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset_Raw"
vcftools --gzvcf $SNPs --keep ${OUTPUTDIR}/SAM_SampleNames288.txt --recode --recode-INFO-all --out ${OUTPUTDIR}/SAM_lines/UBC_Dataset_SAMlines
# After filtering, kept 288 out of 1293 Individuals

```