#!/bin/bash

mitoCN_PATH=$(cd .. && pwd -P)
# mitoCN_PATH=/path/to/mitoCN

### download example data to in to the example folder using
# wget https://zenodo.org/record/7964357/files/data.tar.gz
# mv data.tar.gz ${mitoCN_PATH}/example/data.tar.gz

DATA=${mitoCN_PATH}/example/data.tar.gz
if [ ! -f "$DATA" ]; then
    echo "data.tar.gz does not exist."; 
    echo "please download data.tar.gz into the example folder."; 
    exit 0;
fi

tar -xf $DATA

# CRAM example with hg38 ref 
CRAM=${mitoCN_PATH}/example/data/sample.cram
FASTA=${mitoCN_PATH}/example/data/human38_select.fa
out=${mitoCN_PATH}/example

# estimate mtDNA-CN using chr20
bash ${mitoCN_PATH}/mitoCN.sh -d ${mitoCN_PATH} -f $CRAM -a $FASTA -v hg38 -m chrM -k chr20 -o $out
