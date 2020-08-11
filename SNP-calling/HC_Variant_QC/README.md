# Quality Control checks for high-confidence subset

### Navigation: Jump to Section

- [Compare with Truth Set](#compare-with-truth-set)

--

## Compare with Truth Set



```bash
#To use bcftools need to bgzip and index vcf files
bgzip $HC_Subset
tabix Sunflower_SAM_SNP_Calling_high_confidence_subset.vcf.gz

#define variables
HC_Subset="/scratch/eld72413/SAM_seq/results2/VCF_results_new/Create_HC_Subset/Sunflower_SAM_SNP_Calling_high_confidence_subset.vcf.gz"
Truth_Set="/scratch/eld72413/SNParray/MapUniqueSNP_idt90_rename_rmContigs.vcf.gz"

module load BCFtools/1.10.2-GCC-8.3.0
bcftools stats $Truth_Set $HC_Subset
```
Number of SNPs shared between Truth Set + HC Subset: 1077
Number of SNPs only in HC Subset: 9363697
Number of SNPs only "Truth" Set: 5447
