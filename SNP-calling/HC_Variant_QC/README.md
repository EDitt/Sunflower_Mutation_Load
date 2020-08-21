# Quality Control checks for high-confidence subset

### Navigation: Jump to Section

- [Overview](#overview)
- [Compare with Truth Set](#compare-with-truth-set)
- [Compare with XRQ SNPs](#compare-with-xrq-snps)
--

## Overview

#### Run 1
Filtering parameters:  
	- DP per sample cutoff: 5  
	- GQ cutoff: 20  
	- Max "bad" sites: 58  
	- Max heterozygous sites: 259  
	- QUAL score cutoff: 40

## Jupyter Notebook Exploration

Port Jupyter Notebook into Sapelo

setup-

convert VCF to h5
Used vcf_toh5.py script from sequence handling
```bash
module load Anaconda3/2020.02
pip install scikit-allel[full]
qlogin

python vcf_to_h5.py <vcf>
```

```bash
#interactive job
qlogin

# set up environment
module load Anaconda3/2020.02
pip install scikit-allel[full]

NOTEBOOKPORT=8656 #random value between 8000 and 10000

IPUSED=$(hostname -i)

echo "NOTEBOOKPORT is " $NOTEBOOKPORT
echo "IPUSED is " $IPUSED
jupyter-notebook --port $NOTEBOOKPORT --ip=$IPUSED --no-browser

# establish another ssn tunnel into port NOTEBOOKPORT
ssh -N -L 8656:${IPUSED}:8656 eld72413@sapelo2.gacrc.uga.edu

# http://localhost:8656

```


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

## Compare with XRQ SNPs
##### Previous SNP set (called against XRQ, re-aligned against Ha412v2)

```bash
module load BCFtools/1.10.2-GCC-8.3.0

bgzip remappedHa412HO.vcf
tabix remappedHa412HO.vcf.gz

Re_mappedXRQ=/scratch/eld72413/SNPfiles/easyGWAS_Mar5_XRQremapHa412_Burke/remappedHa412HO.vcf.gz

bcftools stats $HC_Subset $Re_mappedXRQ
```

Number of SNPs shared between XRQ mapped SNPs + HC Subset: 238,985
Number of SNPs only in HC Subset: 9,125,789
Number of SNPs only XRQ mapped SNPs: 1,916,391

