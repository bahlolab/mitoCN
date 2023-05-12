#!/bin/bash

# CRAM example with hg38 ref 
# estimate mtDNA-CN using k500 (select 500 bins from each austosome)
mitoCN_PATH=/stornext/Bioinf/data/lab_bahlo/users/wang.lo/git/mitoCN
CRAM=${mitoCN_PATH}/data/sample1.cram
FASTA=${mitoCN_PATH}/data/Homo_sapiens_assembly38.fasta
out=${mitoCN_PATH}/example

bash ${mitoCN_PATH}/mitoCN.sh -d ${mitoCN_PATH} -f $CRAM -a $FASTA -v hg38 -m chrM -o $out

# BAM example with hg19 ref 
# estimate mtDNA-CN using chromosome 20
mitoCN_PATH=/stornext/Bioinf/data/lab_bahlo/users/wang.lo/git/mitoCN
BAM=${mitoCN_PATH}/data/sample2.bam
out=${mitoCN_PATH}/example

bash ${mitoCN_PATH}/mitoCN.sh -d ${mitoCN_PATH} -f $BAM -v hg19 -m chrM -k chr20 -o $out

#SBATCH --job-name=test
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB

