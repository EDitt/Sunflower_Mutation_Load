
## Get stats about missense vs. synonymous SNPs

#### Frequency distribution

```bash
# list of positions
grep -v "#" fullsam_missense.txt | awk '{print $2}' | awk -F ":" -v OFS="\t" '{$1=$1; print $0}' > fullsam_missense_positions.txt
grep -v "#" fullsam_synon.txt | awk '{print $2}' | awk -F ":" -v OFS="\t" '{$1=$1; print $0}' > fullsam_synon_positions.txt
SAMVCF=/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset/Annuus.tranche90.snp.fullsam.90.bi.remappedHa412HO_reheader.vcf.gz

module load VCFtools/0.1.15-foss-2016b-Perl-5.24.1
vcftools --gzvcf $SAMVCF --freq --positions fullsam_missense_positions.txt --out SAM_MISSENSE
# After filtering, kept 119,332 out of a possible 2155376 Sites
vcftools --gzvcf $SAMVCF --freq --positions fullsam_synon_positions.txt --out SAM_SYNON
# After filtering, kept 184,668 out of a possible 2155376 Sites

grep -v "#" wildenv_missense.txt | awk '{print $2}' | awk -F ":" -v OFS="\t" '{$1=$1; print $0}' > wildenv_missense_positions.txt
grep -v "#" wildenv_synon.txt | awk '{print $2}' | awk -F ":" -v OFS="\t" '{$1=$1; print $0}' > wildenv_synon_positions.txt
WILDANN=/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset/Annuus.tranche90.snp.env.90.bi.remappedHa412HO_reheader.vcf.gz

vcftools --gzvcf $WILDANN --freq --positions wildenv_missense_positions.txt --out WILD_MISSENSE
# After filtering, kept 306,904 out of a possible 4,882,321 Sites
vcftools --gzvcf $WILDANN --freq --positions wildenv_synon_positions.txt --out WILD_SYNON
# After filtering, kept 468,294 out of a possible 4,882,321 Sites

```

#### Subset VCF files

```bash
module load VCFtools/0.1.15-foss-2016b-Perl-5.24.1
module load BCFtools/1.10.2-GCC-8.3.0

### SAM lines
SAMVCF=/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset/Annuus.tranche90.snp.fullsam.90.bi.remappedHa412HO_reheader.vcf.gz
SAM_MISSENSE_POS=/scratch/eld72413/NSFproj/VEP/Missense/fullsam_missense_positions.txt

vcftools --gzvcf $SAMVCF --positions $SAM_MISSENSE_POS --recode --recode-INFO-all --out SAM_MISSENSE
bgzip SAM_MISSENSE.recode.vcf

SAM_SYNON_POS=/scratch/eld72413/NSFproj/VEP/Synon/fullsam_synon_positions.txt

vcftools --gzvcf $SAMVCF --positions $SAM_SYNON_POS --recode --recode-INFO-all --out SAM_SYNON
bgzip SAM_SYNON.recode.vcf


### Wild lines

WILDANN=/scratch/eld72413/NSFproj/PublishedSNPs/UBC_Dataset/Annuus.tranche90.snp.env.90.bi.remappedHa412HO_reheader.vcf.gz



```

#### Run through VeP to get figures showing density per chromosome for different variant classes