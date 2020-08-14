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

# used tmux session for longer download:
aws s3 cp s3://ubc-sunflower-genome/cohorts/ha412v2/wgs_all/samples.json .
```