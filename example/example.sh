#!/bin/bash

# CRAM example with hg38 ref 
# estimate mtDNA-CN using chr20
mitoCN_PATH=/stornext/Bioinf/data/lab_bahlo/users/wang.lo/git/mitoCN
tar -xf ${mitoCN_PATH}/example/data.gz
CRAM=${mitoCN_PATH}/example/data/sample.cram
FASTA=${mitoCN_PATH}/example/data/human38_select.fa
out=${mitoCN_PATH}/example

bash ${mitoCN_PATH}/mitoCN.sh -d ${mitoCN_PATH} -f $CRAM -a $FASTA -v hg38 -m chrM -k chr20 -o $out
