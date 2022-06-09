# Numbers in each variant class

## Navigation: Jump to Section

- [VeP annotation classes represented](#)
- [0. Parse deleterious predictions](#0.-parse-deleterious-predictions)

## VeP annotation classes represented:

```bash
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

# VeP categories:
grep -v "#" $VEP | awk '{print $7}' | sort -u

```

In order of severity (according to VeP):
1.) Splice donor/acceptor variant
		splice_acceptor_variant
		splice_donor_variant
2.) Stop/Start Lost/Gained
		start_lost
		start_lost,splice_region_variant
		stop_gained
		stop_gained,splice_region_variant
		stop_lost
		stop_lost,splice_region_variant
3.) Missense Variant
		missense_variant
4.) Synonymous (includes splice region variant and stop retained variants)
		missense_variant,splice_region_variant (because this wasn't in the group tested by BAD_Mutations)
		stop_retained_variant
		synonymous_variant
		splice_region_variant,stop_retained_variant
		splice_region_variant,synonymous_variant
		splice_region_variant,3_prime_UTR_variant
		splice_region_variant,5_prime_UTR_variant
		splice_region_variant,intron_variant
5.) Non-coding (includes 5' UTR, 3' UTR, intron, upstream, downstream, intergenic)
		3_prime_UTR_variant
		5_prime_UTR_variant
		downstream_gene_variant
		intergenic_variant
		intron_variant
		upstream_gene_variant


## 0.) Parse deleterious predictions
- The number of deleterious and tolerated sites were identified and duplicates excluded from variant lists (see `2.BAD_Mutations/Post_processing.md`)

### Deleterious
```bash
srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
source /home/eld72413/DelMut/Sunflower_Mutation_Load/BAD_Mutations/Variant_analyses/config.sh

grep "Tolerated" ${DSNP_DATA} | wc -l # 557,317
grep "Deleterious" ${DSNP_DATA} | wc -l # 87,891

# how many duplicates?
awk 'NR>1 {print $2}' $DSNP_DATA | wc -l # 645,208
awk 'NR>1 {print $2}' $DSNP_DATA | sort -u | wc -l # 641,514

awk -F'\t' '{print NF; exit}' ${DSNP_DATA} # 39 columns

## deleterious positions in reference:
awk '{if ($39=="Reference_deleterious") {print $0}}' ${DSNP_DATA} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' | wc -l # 11,796
awk '{if ($39=="Reference_deleterious") {print $0}}' ${DSNP_DATA} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' > $OUT_DIR/IntermediateFiles/Reference_DelPositions.txt # 11,796 no duplicates


# alternate deleterious:
awk '{if ($39=="Alternate_deleterious") {print $0}}' ${DSNP_DATA} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' | wc -l # 76,095

awk '{if ($39=="Alternate_deleterious") {print $0}}' ${DSNP_DATA} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' | sort -u > $OUT_DIR/IntermediateFiles/Alternate_DelPositions.txt # 76,016

# Total dSNPs
cat $OUT_DIR/IntermediateFiles/Reference_DelPositions.txt $OUT_DIR/IntermediateFiles/Alternate_DelPositions.txt | wc -l # 87,812

cat $OUT_DIR/IntermediateFiles/Reference_DelPositions.txt $OUT_DIR/IntermediateFiles/Alternate_DelPositions.txt |sort -u > $OUT_DIR/SupportingFiles/FinalPositionFiles/AllDel_Positions.txt 
# 87794 (87812 - 18)

```

### Tolerated
```bash
awk '{if ($39=="Tolerated") {print $0}}' ${DSNP_DATA} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' | wc -l # 556,577 (740 duplicate positions represented in deleterious set already removed - see 2.BAD_Mutations/Post_processing.md)

awk '{if ($39=="Tolerated") {print $0}}' ${DSNP_DATA} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' | sort -u > $OUT_DIR/IntermediateFiles/ToleratedPositions.txt # 553,720

# check
comm -12 --check-order $OUT_DIR/IntermediateFiles/ToleratedPositions.txt $OUT_DIR/SupportingFiles/FinalPositionFiles/AllDel_Positions.txt | wc -l # 0

cp $OUT_DIR/IntermediateFiles/ToleratedPositions.txt $OUT_DIR/SupportingFiles/FinalPositionFiles
```

### Total number that were aligned and tested
```bash
# number that were aligned + tested
# combine with tolerated
cat $OUT_DIR/SupportingFiles/FinalPositionFiles/AllDel_Positions.txt $OUT_DIR/IntermediateFiles/ToleratedPositions.txt | sort -u  > $OUT_DIR/IntermediateFiles/Missense_Tested.txt # 641,514

```
After removing duplicate sites:
- Total Deleterious: 87,794 (87,812 - 18)
	- Reference allele deleterious: 11,796
	- Alternate allele deleterious: 76,016
- Total Tolerated: 553,720

- Total tested: 641,514

Note: there are 18 positions that are the same in the alternate deleterious & reference deleterious sets (not removed)

Duplicates removed:
- 79 duplicates within alternate deleterious category
- 740 duplicates that were represented in deleterious & tolerated categories (removed from tolerated)
- 2,857 duplicates within tolerated category

## 1.) Splice donor/acceptor variant

```bash
# 1.)
grep -v "#" $VEP | awk '{if ($7=="splice_acceptor_variant" || $7=="splice_donor_variant") {print $2}}' | awk -F":" '{print $1"\t"$2}' | sort -u > $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt # 6907

# remove tested missense SNPs
comm -23 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/Missense_Tested.txt | wc -l # 6851

comm -23 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/Missense_Tested.txt > $OUT_DIR/SupportingFiles/FinalPositionFiles/Splice_AcceptorDonor_Nodups.txt

```

## 2.) Start/Stop Lost/Gained

```bash
# 2.)
grep -v "#" $VEP | awk '{if ($7~"lost" || $7~"gained") {print $2}}' | awk -F":" '{print $1"\t"$2}' | sort -u > $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt # 23,734

### VeP 'most severe' number is: 23,732 (so most likely 2 are also annotated as splice acceptor/donor variant)
comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt | wc -l # 2 - yes

comm -13 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt > $OUT_DIR/IntermediateFiles/StartStop_LostGained_noDups.txt # 23,732

# remove tested missense mutations
comm -23 --check-order $OUT_DIR/IntermediateFiles/StartStop_LostGained_noDups.txt $OUT_DIR/IntermediateFiles/Missense_Tested.txt > $OUT_DIR/SupportingFiles/FinalPositionFiles/StartStop_LostGained_noDups.txt # 23,383
```

concatenate the two classes & with tested missense mutations
```bash
cat $OUT_DIR/SupportingFiles/FinalPositionFiles/Splice_AcceptorDonor_Nodups.txt $OUT_DIR/SupportingFiles/FinalPositionFiles/StartStop_LostGained_noDups.txt | sort -u > $OUT_DIR/IntermediateFiles/severe12.txt # 30,234

cat $OUT_DIR/IntermediateFiles/severe12.txt $OUT_DIR/IntermediateFiles/Missense_Tested.txt | sort -u > $OUT_DIR/IntermediateFiles/severe12_missenseTested.txt # 671,748

```

## 3.) Missense Variants

```bash
# missed due to double annotation
grep -v "#" $VEP | awk '{if ($7=="missense_variant,splice_region_variant") {print $2}}' | awk '{$1=$1}1' FS=':' OFS='\t' | sort -u > $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated.txt # 8555

# all missense
grep -v "#" $VEP | awk '{if ($7=="missense_variant") {print $2}}' | awk '{$1=$1}1' FS=':' OFS='\t' | sort -u > $OUT_DIR/IntermediateFiles/All_missense.txt # 699,805

# any shared between these two annotations?
comm -12 --check-order $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated.txt $OUT_DIR/IntermediateFiles/All_missense.txt | wc -l # 49

# list of "missense, splice region variant" SNPs that are not in "All missense"
comm -23 --check-order $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated.txt $OUT_DIR/IntermediateFiles/All_missense.txt > $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated_NoDups.txt # 8506

# in more severe class?
comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated_NoDups.txt | wc -l # 11
comm -12 --check-order $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated_NoDups.txt | wc -l # 0

# which missense snps are in more severe class
comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/All_missense.txt | wc -l # 62
comm -12 --check-order $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt $OUT_DIR/IntermediateFiles/All_missense.txt | wc -l # 379

# remove all previous classes
comm -13 --check-order $OUT_DIR/IntermediateFiles/severe12_missenseTested.txt $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated_NoDups.txt > $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated_NoDups2.txt

```
there were 8495 missense SNPs that were not tested due to being annotated as "missense_variant,splice_region_variant" & not in a more severe category
(49 annotated as "missense" in a duplicate position, 11 annotated as splice acceptor/donor)


### Un-alignable - some missense variants were not able to be aligned and thus have no tolerated/deleterious annotation
```bash

# number that were not aligned and tested
comm -23 --check-order $OUT_DIR/IntermediateFiles/All_missense.txt $OUT_DIR/IntermediateFiles/Missense_Tested.txt | sort -u > $OUT_DIR/IntermediateFiles/Missense_UnAligned.txt #58,291

comm -12 --check-order $OUT_DIR/SupportingFiles/FinalPositionFiles/Splice_AcceptorDonor_Nodups.txt $OUT_DIR/IntermediateFiles/Missense_UnAligned.txt | wc -l # 6

comm -12 --check-order $OUT_DIR/SupportingFiles/FinalPositionFiles/StartStop_LostGained_noDups.txt $OUT_DIR/IntermediateFiles/Missense_UnAligned.txt | wc -l # 30

comm -13 --check-order $OUT_DIR/IntermediateFiles/severe12_missenseTested.txt $OUT_DIR/IntermediateFiles/Missense_UnAligned.txt > $OUT_DIR/IntermediateFiles/Missense_UnAligned_noDups.txt # 58,255
```

Combine the two groups of 'missed' missense mutations
```bash
cat $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated_NoDups2.txt $OUT_DIR/IntermediateFiles/Missense_UnAligned_noDups.txt | sort -u > $OUT_DIR/SupportingFiles/FinalPositionFiles/Missense_Other.txt # 66,750

# combine all severe variant classes
cat $OUT_DIR/IntermediateFiles/severe12_missenseTested.txt $OUT_DIR/SupportingFiles/FinalPositionFiles/Missense_Other.txt | sort -u > $OUT_DIR/IntermediateFiles/severe12_missense.txt # 738,498
```

## 4.) Synonymous Variants
```bash
grep -v "#" $VEP | awk '{if ($7~"splice_region_variant" || $7~"synonymous_variant" || $7~"stop_retained_variant") {print $2}}' | awk -F":" '{print $1"\t"$2}' | sort -u > $OUT_DIR/IntermediateFiles/Synonymous_positions.txt # 943,011

# remove those represented in a more severe variant class
comm -13 --check-order $OUT_DIR/IntermediateFiles/severe12_missense.txt $OUT_DIR/IntermediateFiles/Synonymous_positions.txt > $OUT_DIR/SupportingFiles/FinalPositionFiles/Synonymous_Nodups.txt # 927,677


# combine all previous variant classes
cat $OUT_DIR/IntermediateFiles/severe12_missense.txt $OUT_DIR/SupportingFiles/FinalPositionFiles/Synonymous_Nodups.txt | sort -u > $OUT_DIR/IntermediateFiles/All_Coding.txt # 1,666,175
```

## 5.) Non-coding Variants

```bash
grep -v "#" $VEP | awk '{if ($7~"UTR" || $7~"intergenic" || $7~"downstream" || $7~"upstream" || $7~"intron") {print $2}}' | awk -F":" '{print $1"\t"$2}' | sort -u > $OUT_DIR/IntermediateFiles/NonCoding.txt # 36,336,841

# remove those represented in a more severe variant class
comm -13 --check-order $OUT_DIR/IntermediateFiles/All_Coding.txt $OUT_DIR/IntermediateFiles/NonCoding.txt > $OUT_DIR/SupportingFiles/FinalPositionFiles/NonCoding_Nodups.txt # 35,453,937
```