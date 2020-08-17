Download SNP set from UBC Sunflower Genome aws bucket

install modules needed
```bash
module load Anaconda3/2020.02
pip install --upgrade pip
pip install boto3 awscli

# set path environmental variable
export PATH=$PATH:/home/eld72413/.local/bin
```

configure AWS
```bash
aws configure
# put in access key ID, secret access key, us-west-2 as region, json as default output format
```

download files
```bash
aws s3 cp s3://ubc-sunflower-genome/cohorts/ha412v2/wgs_all/sample-names.tsv .

#list of files to download
aws s3 ls s3://ubc-sunflower-genome/haploblocks/processed_snps/
aws s3 ls s3://ubc-sunflower-genome/haploblocks/processed_snps/all_ann1/

aws s3 cp s3://ubc-sunflower-genome/haploblocks/processed_snps/listing.sorted.txt .

aws s3 ls s3://ubc-sunflower-genome/haploblocks/processed_snps/all_ann1/fullsam/
aws s3 cp s3://ubc-sunflower-genome/haploblocks/processed_snps/all_ann1/fullsam/fullsam_sample.list . #the 287 samples
aws s3 cp s3://ubc-sunflower-genome/haploblocks/processed_snps/all_ann1/fullsam/sam_sample.list . #264

aws s3 cp s3://ubc-sunflower-genome/haploblocks/processed_snps/all_ann1/sample_sets/env_sample.list . #719 samples

aws s3 cp s3://ubc-sunflower-genome/haploblocks/processed_snps/all_ann1/sample_sets/herb_sample.list . #323 samples

aws s3 cp s3://ubc-sunflower-genome/haploblocks/processed_snps/all_ann1/sample_sets/gwas_sample.list . #614 samples

aws s3 cp s3://ubc-sunflower-genome/haploblocks/processed_snps/all_ann1/sample_sets/ww_samples.list . #163 samples

https://ubc-sunflower-genome.s3-us-west-2.amazonaws.com/haploblocks/processed_snps/listing.sorted.txt  
# used tmux session for longer download:
aws s3 cp s3://ubc-sunflower-genome/cohorts/ha412v2/wgs_all/samples.json .
```

### Wild Annuus data from Todesco et al.
- 719 re-sequenced individuals (listed in "Coverage and analyses")
	- 265 were "samples from other studies"