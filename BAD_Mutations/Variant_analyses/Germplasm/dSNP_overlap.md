# Shared vs. private dSNPs for the different heterotic groups

### Script
# test script
```bash
sbatch --export=group='',alt_vcf='',ref_vcf='',OUTPUT_DIR'='' VCF_SubsetByGenotype.sh
```

### Setup
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
```

### Groupings
```bash
# test code in script with RHA oil
group="RHA-Oil"
alt_vcf=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Altdeleterious.vcf
ref_vcf=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/SAM_Refdeleterious.vcf
OUTPUT_DIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups

awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l #66

# Number of alt alleles: # 28,487
# Number of ref alleles: # 11,694

# concatenated file:
grep -v "#" RHA-Oil_AllDel.vcf | wc -l # 40,181
```

# compare number same/different
```bash
OUTPUTDIR=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/intersect_tests
HA_oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/HA-Oil_AllDel.vcf
RHA_oil=/scratch/eld72413/SAM_seq/BAD_Mut_Files/Results/Groups/RHA-Oil_AllDel.vcf

bgzip -c ${HA_oil} > ${HA_oil}.gz
tabix -p vcf ${HA_oil}.gz
bgzip -c ${RHA_oil} > ${RHA_oil}.gz
tabix -p vcf ${RHA_oil}.gz

### needs to be compressed!
# -c all means that it considers positions identical regarless of whether ALT allele matches or not
# should not matter here because all sites biallelic
bcftools isec -p ${OUTPUTDIR}/Oil ${HA_oil}.gz ${RHA_oil}.gz

grep -v "#" 0000.vcf | wc -l # 10,505 (private to HA-oil)
grep -v "#" 0001.vcf | wc -l # 10,599 (private to RHA-oil)
grep -v "#" 0002.vcf | wc -l # 29,582 HA-oil shared by both
grep -v "#" 0003.vcf | wc -l # 29,582 RHA-oil shared by both

# totals 50,686
# HA Oil: 40,087 (agrees with these numbers)
# RHA Oil: 40,181 (agrees with these numbers)
# different with the -c all flag?
bcftools isec -p ${OUTPUTDIR}/Oil2 -c all ${HA_oil}.gz ${RHA_oil}.gz

### numbers exactly the same
```

#### HA-oil
```bash
# start with HA-oil:

SAM_info=/home/eld72413/DelMut/Sunflower_Mutation_Load/SNP-calling/All_SAM_Info.csv

# last two lines are not included in vcf
group="HA-Oil"

# need comma separated list of samples
# need to remove the 2 samples that are not in the VCF
genotypes=$(awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | paste -sd,)
awk -v var="$group" -F',' '{if ($8==var && $9!="HA412" && $9!="NA") {print $9}}' $SAM_info | wc -l # 63

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

bcftools concat -a ${OUTPUT_DIR}/${group}_AltDel.vcf ${OUTPUT_DIR}/${group}_RefDel2.vcf -o ${OUTPUT_DIR}/${group}_AllDel.vcf
# to use -a flag, files need to be gzipped?

bgzip -c ${OUTPUT_DIR}/${group}_AltDel.vcf > ${OUTPUT_DIR}/${group}_AltDel.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${group}_AltDel.vcf.gz
bgzip -c ${OUTPUT_DIR}/${group}_RefDel2.vcf > ${OUTPUT_DIR}/${group}_RefDel2.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${group}_RefDel2.vcf.gz

# also needs index
bcftools concat -a ${OUTPUT_DIR}/${group}_AltDel.vcf.gz ${OUTPUT_DIR}/${group}_RefDel2.vcf.gz -Ov -o ${OUTPUT_DIR}/${group}_AllDel.vcf
grep -v "#" HA-Oil_AllDel.vcf | wc -l # 40,087 <- this gives the right number!

# also test with MergeVcfs (Picard)
module load picard/2.21.6-Java-11
java -jar /apps/eb/picard/2.21.6-Java-11/picard.jar MergeVcfs \
          I=${OUTPUT_DIR}/${group}_AltDel.vcf \
          I=${OUTPUT_DIR}/${group}_RefDel2.vcf \
          O=output_variants.vcf
grep -v "#" output_variants.vcf | wc -l #0

# count variants in gzipped vcf?
bcftools view -Ov ${OUTPUT_DIR}/${group}_AltDel.vcf.gz | grep -v "#" | wc -l # 28,311

#### tests
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




```