# Shared vs. private dSNPs for the different heterotic groups

### Groupings

```bash
# start with HA-oil:

SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv

# last two lines are not included in vcf
group="HA-Oil"

# need comma separated list of samples
# need to remove the 2 samples that are not in the VCF
genotypes=$(awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | paste -sd,)
awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l

# for dSNPs that are the alternate allele
alt_vcf=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf
grep -v "#" $alt_vcf | wc -l # 76,016

OUTPUT_DIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups

module load BCFtools/1.10.2-GCC-8.3.0

bcftools query -l ${alt_vcf}

bcftools view -Ou --samples ${genotypes} ${alt_vcf} | bcftools filter -i 'COUNT(GT="alt") > 0' -o ${OUTPUT_DIR}/${group}_AltDel.vcf

grep -v "#" HA-Oil_AltDel.vcf | wc -l #28,311

# does this number differ if I use GT="AA" (require homozygous)?

bcftools view -Ou --samples ${genotypes} ${alt_vcf} | bcftools filter -i 'COUNT(GT="AA") > 0' -o ${OUTPUT_DIR}/${group}_AltDel_hom.vcf

grep -v "#" HA-Oil_AltDel_hom.vcf | wc -l #15,979

bcftools view -Ou --samples ${genotypes} ${alt_vcf} | bcftools filter -i 'COUNT(GT="AA" | GT="het") > 0' -o ${OUTPUT_DIR}/${group}_AltDel2.vcf

grep -v "#" HA-Oil_AltDel2.vcf | wc -l # 28,311


# for dSNPs that are the reference allele
ref_vcf=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf
grep -v "#" $ref_vcf | wc -l # 11,796

bcftools view -Ou --samples ${genotypes} ${ref_vcf} | bcftools filter -i 'COUNT(GT="ref") > 0' -o ${OUTPUT_DIR}/${group}_RefDel.vcf
grep -v "#" HA-Oil_RefDel.vcf | wc -l # 11,774

bcftools view -Ou --samples ${genotypes} ${ref_vcf} | bcftools filter -i 'COUNT(GT="RR") > 0' -o ${OUTPUT_DIR}/${group}_RefDel_hom.vcf
grep -v "#" HA-Oil_RefDel_hom.vcf | wc -l # 11,774

bcftools view -Ou --samples ${genotypes} ${ref_vcf} | bcftools filter -i 'COUNT(GT="ref"| GT="het") > 0' -o ${OUTPUT_DIR}/${group}_RefDel2.vcf
grep -v "#" HA-Oil_RefDel2.vcf | wc -l # 11,776 - need this for the reference alleles (if I'm allowing hets)


### combine the files:
bcftools concat ${OUTPUT_DIR}/${group}_AltDel.vcf ${OUTPUT_DIR}/${group}_RefDel2.vcf -o ${OUTPUT_DIR}/${group}_AllDel.vcf
# Checking the headers and starting positions of 2 files
# Concatenating /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/HA-Oil_AltDel.vcf	0.315904 seconds
# Concatenating /scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/HA-Oil_RefDel2.vcf
#The chromosome block Ha412HOChr00c00004 is not contiguous, consider running with -a.

grep -v "#" HA-Oil_AllDel.vcf | wc -l # 28,080 (expecting 40,087, difference of 12,007)


bcftools concat $alt_vcf $ref_vcf -o ${OUTPUT_DIR}/AllDel_TEST.vcf
bcftools concat $alt_vcf $alt_vcf -o ${OUTPUT_DIR}/AllDel_TEST2.vcf
grep -v "#" AllDel_TEST.vcf | wc -l # 75,982 (should be 87,812) (difference of 11,830)
grep -v "#" AllDel_TEST2.vcf | wc -l # 75,982
cat $alt_vcf $ref_vcf > ${OUTPUT_DIR}/AllDel_TEST3.vcf
grep -v "#" AllDel_TEST3.vcf | wc -l # 87,812
cat $alt_vcf $alt_vcf > ${OUTPUT_DIR}/AllDel_TEST4.vcf
grep -v "#" AllDel_TEST4.vcf | wc -l  # 152,032

### why is the concatenated vcf for HA oil not the same as the sum of the 2?

bcftools query -l ${OUTPUT_DIR}/${group}_AltDel.vcf | wc -l # 63
bcftools query -l ${OUTPUT_DIR}/${group}_RefDel2.vcf | wc -l # 63


##############
# for all HAs
awk -F',' '{if ($8~"HA" && $8!~"RHA") {print $9}}' $SAM_info

# for all RHAs
awk -F',' '{if ($8~"RHA") {print $9}}' $SAM_info





bcftools filter -i 'AC > 0'  # what is the difference?


#GT = "AA" if I want homozygous only





bcftools filter -i 'COUNT(GT="het")/(N_SAMPLES-N_MISSING) < 0.2' ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_missing_filtered.vcf -o ${OUTPUT_DIR}/Intermediates/${OUT_PREFIX}_HETFiltered.vcf
```