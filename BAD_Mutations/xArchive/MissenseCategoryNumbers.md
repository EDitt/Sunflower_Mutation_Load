# Missense Numbers

## Navigation: Jump to Section

- [Summay of missense SNPs](#summary-of-missense-snps-tested)
- [Parse deleterious predictions](#parse-deleterious-predictions)

## Summary of Missense SNPs tested

For Pie graph will use the numbers from VeP "most severe consequence", however need to further break down missense category:

### Summary of results
- 707,859 missense mutations (as 'most severe' consequence)
	- 699,805 missense mutations were tested with BAD_Mutations:
		- 553,702 were annotated as tolerated (353 have more severe consequence)
		- 87,794 were annotated as deleterious (52 have more severe consequence)
		- 58,291 were unable to be aligned (36 have more severe consequence)
	- 8506 were missed due to having annotation as "missense_variatn, splice_region_variant" (11 have more severe consequence)

```bash
wc -l All_SNP_Info.txt # 36,708,693
```
### From VeP:

```bash
# VeP categories:
grep -v "#" $VEP | awk '{print $7}' | sort -u

### other categories with 'missense': 
# missense_variant,splice_region_variant

grep -v "#" $VEP | awk '{print $2}' | wc -l # 43,019,114
grep -v "#" $VEP | awk '{print $2}' | sort -u | wc -l # 37,120,112 (same as VeP "most severe" total)

```

```bash
# number tested with BAD_Mutations (used =='missense' in python convert.py script)
grep -v "#" $VEP | awk '{if ($7=="missense_variant") {print $2}}' | wc -l # 704,075
grep -v "#" $VEP | awk '{if ($7=="missense_variant") {print $2}}' | sort -u | wc -l # 699,805 (4,270 duplicates)

# all missense:
grep -v "#" $VEP | awk '{if ($7~"missense_variant") {print $2}}' | wc -l # 712,631
grep -v "#" $VEP | awk '{if ($7~"missense_variant") {print $2}}' | sort -u | wc -l # 708,311

# missing from BAD_Mutations list of tested:
grep -v "#" $VEP | awk '{if ($7=="missense_variant,splice_region_variant") {print $2}}' | wc -l # 8556

```

8556 were not tested due to having an annotation of: "missense_variant,splice_region_variant"


### The number of predicted consequences:

```bash
wc -l SAM_SNP_BadMut_Summary # 704,075 predictions tested
awk '{print $4}' SAM_SNP_BadMut_Summary | sort -u | wc -l # 699,805
```
The number of sites tested with BAD_Mutations was: 704,075 (unique SNPs: 699,805)


### Predictions Returned:

```bash
awk 'NR>1 {print $38}' dsnp_data_Polarized.table | wc -l # 645208 predictions returned
```
### Numbers of dSNPs minus number untested?

### More severe consequences:
According to VeP, these categories have more severe consequences and are not counted in missense number for "most severe":

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

### Lists of more severe consequences:
```bash
# 1.)
grep -v "#" $VEP | awk '{if ($7=="splice_acceptor_variant" || $7=="splice_donor_variant") {print $2}}' | sort -u > $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt # 6907

# 2.)
grep -v "#" $VEP | awk '{if ($7~"lost" || $7~"gained") {print $2}}' | sort -u > $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt # 23,734
### VeP 'most severe' number is: 23,732 (so most likely 2 are also annotated as splice acceptor/donor variant)

comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt | wc -l # 2 - yes
```

### Lists of missense categories
```bash
# missed due to double annotation
grep -v "#" $VEP | awk '{if ($7=="missense_variant,splice_region_variant") {print $2}}' | sort -u > $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated.txt # 8555

# all missense
grep -v "#" $VEP | awk '{if ($7=="missense_variant") {print $2}}' | sort -u > $OUT_DIR/IntermediateFiles/All_missense.txt # 699,805

# any shared between these two annotations?
comm -12 --check-order $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated.txt $OUT_DIR/IntermediateFiles/All_missense.txt | wc -l # 49

# list of "missense, splice region variant" SNPs that are not in "All missense"
comm -23 --check-order $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated.txt $OUT_DIR/IntermediateFiles/All_missense.txt > $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated_NoDups.txt # 8506


comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated_NoDups.txt | wc -l # 11
comm -12 --check-order $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt $OUT_DIR/IntermediateFiles/Missense_DoubleAnnotated_NoDups.txt | wc -l # 0


comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/All_missense.txt | wc -l # 62
comm -12 --check-order $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt $OUT_DIR/IntermediateFiles/All_missense.txt | wc -l # 379


```
there were 8495 missense SNPs that were not tested due to being annotated as "missense_variant,splice_region_variant" & not in a more severe category



### Numbers in various allele classes (excluding duplicates)
- The number of deleterious and tolerated sites were identified and duplicates excluded from variant lists (see `2.BAD_Mutations/Post_processing.md`)


# Parse deleterious predictions

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
awk '{if ($39=="Reference_deleterious") {print $0}}' ${DSNP_DATA} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' > $OUT_DIR/IntermediateFiles/Reference_DelPositions.txt
wc -l $OUT_DIR/IntermediateFiles/Reference_DelPositions.txt # 11,796
awk '{print $0}' $OUT_DIR/IntermediateFiles/Reference_DelPositions.txt | sort -u | wc -l # 11,796 no duplicates

mv $OUT_DIR/IntermediateFiles/Reference_DelPositions.txt /scratch/eld72413/SAM_seq/dSNP_results/SupportingFiles/FinalPositionFiles

# alternate deleterious:
awk '{if ($39=="Alternate_deleterious") {print $0}}' ${DSNP_DATA} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' > $OUT_DIR/IntermediateFiles/Alternate_DelPositions.txt
wc -l $OUT_DIR/IntermediateFiles/Alternate_DelPositions.txt # 76,095
awk '{print $0}' $OUT_DIR/IntermediateFiles/Alternate_DelPositions.txt | sort -u | wc -l # 76,016
awk '{print $0}' $OUT_DIR/IntermediateFiles/Alternate_DelPositions.txt | sort -u  > $OUT_DIR/SupportingFiles/FinalPositionFiles/Alternate_DelPositionsNoDups.txt
wc -l $OUT_DIR/SupportingFiles/FinalPositionFiles/Alternate_DelPositionsNoDups.txt # 76,016


# Tolerated:
awk '{if ($39=="Tolerated") {print $0}}' ${DSNP_DATA} | awk '{print $25}' | awk '{$1=$1}1' FS=':' OFS='\t' > $OUT_DIR/IntermediateFiles/ToleratedPositions.txt
wc -l $OUT_DIR/IntermediateFiles/ToleratedPositions.txt # 556,577 (740 duplicate positions represented in deleterious set already removed - see 2.BAD_Mutations/Post_processing.md)
awk '{print $0}' $OUT_DIR/IntermediateFiles/ToleratedPositions.txt | sort -u | wc -l # 553,720
awk '{print $0}' $OUT_DIR/IntermediateFiles/ToleratedPositions.txt | sort -u > $OUT_DIR/SupportingFiles/FinalPositionFiles/ToleratedPositionsNoDups.txt
wc -l $OUT_DIR/SupportingFiles/FinalPositionFiles/ToleratedPositionsNoDups.txt # 553,720

```

After removing duplicate sites:
- Total Deleterious: 87,812
	- Reference allele deleterious: 11,796
	- Alternate allele deleterious: 76,016
- Total Tolerated: 553,720

Note: there are 18 positions that are the same in the alternate deleterious & reference deleterious sets (not removed)

Duplicates removed:
- 79 duplicates within alternate deleterious category
- 740 duplicates that were represented in deleterious & tolerated categories (removed from tolerated)
- 2,857 duplicates within tolerated category


### Un-alignable - some missense variants were not able to be aligned and thus have no tolerated/deleterious annotation
```bash
# number that were aligned + tested

cat $OUT_DIR/SupportingFiles/FinalPositionFiles/Alternate_DelPositionsNoDups.txt $OUT_DIR/SupportingFiles/FinalPositionFiles/Reference_DelPositions.txt | sort -u > $OUT_DIR/SupportingFiles/FinalPositionFiles/AllDel_Positions.txt 
# 87794 (87812 - 18)

# combine with tolerated
cat $OUT_DIR/SupportingFiles/FinalPositionFiles/AllDel_Positions.txt $OUT_DIR/SupportingFiles/FinalPositionFiles/ToleratedPositionsNoDups.txt | sort -u | awk '{print $1":"$2}' > $OUT_DIR/IntermediateFiles/Missense_Tested.txt # 641,514

# number that were not aligned and tested

comm -23 --check-order $OUT_DIR/IntermediateFiles/All_missense.txt $OUT_DIR/IntermediateFiles/Missense_Tested.txt | sort -u > $OUT_DIR/IntermediateFiles/Missense_UnAligned.txt
# 58,291

comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/Missense_UnAligned.txt | wc -l # 6
comm -12 --check-order $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt $OUT_DIR/IntermediateFiles/Missense_UnAligned.txt | wc -l # 30

```

More severe consequence shared with deleterious/tolerated
```bash
comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/Missense_Tested.txt | wc -l # 56
comm -12 --check-order $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt $OUT_DIR/IntermediateFiles/Missense_Tested.txt | wc -l # 349

awk '{print $1":"$2}' $OUT_DIR/SupportingFiles/FinalPositionFiles/AllDel_Positions.txt > $OUT_DIR/IntermediateFiles/AllDel_Positions.txt
comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/AllDel_Positions.txt | wc -l # 4
comm -12 --check-order $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt $OUT_DIR/IntermediateFiles/AllDel_Positions.txt | wc -l # 48

awk '{print $1":"$2}' $OUT_DIR/SupportingFiles/FinalPositionFiles/ToleratedPositionsNoDups.txt > $OUT_DIR/IntermediateFiles/ToleratedPositionsNoDups.txt
comm -12 --check-order $OUT_DIR/IntermediateFiles/Splice_AcceptorDonor.txt $OUT_DIR/IntermediateFiles/ToleratedPositionsNoDups.txt | wc -l # 52
comm -12 --check-order $OUT_DIR/IntermediateFiles/StartStop_LostGained.txt $OUT_DIR/IntermediateFiles/ToleratedPositionsNoDups.txt | wc -l # 301
```

```bash
grep -v "#" $VEP | awk '{if ($7!="missense_variant") {print $2}}' | sort -u > $OUT_DIR/IntermediateFiles/No_missense.txt # 36,778,006

```



