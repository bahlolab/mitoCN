#!/bin/bash

# BAM example with hg19 ref 
# estimate mtDNA-CN using chromosome 20
mitoCN_PATH=/stornext/Bioinf/data/lab_bahlo/users/wang.lo/git/mitoCN
BAM=${mitoCN_PATH}/data/sample2.bam
out=${mitoCN_PATH}/example

bash ${mitoCN_PATH}/mitoCN.sh -d ${mitoCN_PATH} -f $BAM -v hg19 -m chrM -k chr20 -o $out