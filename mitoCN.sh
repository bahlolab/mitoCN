#!/bin/bash
#mitoCN 0.1.0

set -euxo pipefail

# mosdepth=/path/to/mosdepth/mosdepth_version
mosdepth=/wehisan/bioinf/lab_bahlo/software/apps/mosdepth/mosdepth_v0.2.9
# module load R

if [ ! -f "$mosdepth" ]; then
    echo "$mosdepth does not exist."; 
    exit 0;
fi

while getopts d:k:m:v:f:a:o: flag
do
    case "${flag}" in
        d) work_dir=${OPTARG};;
        k) region=${OPTARG};;
        m) mt=${OPTARG};;
        v) ref_ver=${OPTARG};;
        f) BAM_file=${OPTARG};;
        a) ref_fasta=${OPTARG};;
        o) out_dir=${OPTARG};;
    esac
done

echo "mitoCN diretory: $work_dir";

if [[ "$region" == "" ]] 
then
    region=k500;
fi

echo "Region: $region";

if [[ "$mt" != "MT" && "$mt" != "chrM" ]]
then
    echo "error: -m must be MT or chrM.";
    exit 1;
else
    if [[ "$mt" == "MT" ]]
    then
        echo "Header: 1,2,...,22,X,Y,$mt";
        X=X; Y=Y;
    else
        echo "Header: chr1,chr2,...,chr22,chrX,chrY,$mt";
        X=chrX; Y=chrY;
    fi
fi
echo "Reference version: $ref_ver";
echo "BAM/CRAM file: $BAM_file";

input_fmt=$(echo "$BAM_file" | awk '{ n=split($BAM_file, arr, "."); print arr[n] }')
echo "Input format: $input_fmt";
echo "fasta file: $ref_fasta";

if [[ "$input_fmt" == "cram" && "$ref_fasta" == "" ]]
then
    echo "error: fasta file must be provided for cram file.";
    exit 2;
fi

ID=$(echo "$BAM_file" | awk '{ n=split($BAM_file, arr, "/"); print arr[n] }' | cut -d '.' -f 1)
echo "Sample ID: $ID";
echo "Output directory: $out_dir";

if [ ! -d "${out_dir}/${ID}" ] 
then
    mkdir -p ${out_dir}/${ID}
else
    echo "Directory ${out_dir}/${ID} exists."
fi

### BED files
BED_file_at=${work_dir}/${ref_ver}/${mt}/${region}.bed
GC_file_at=${work_dir}/${ref_ver}/${mt}/${region}.gc.bed
BED_file_mt=${work_dir}/${ref_ver}/${mt}/chrM.bed
GC_file_mt=${work_dir}/${ref_ver}/${mt}/chrM.gc.bed
BED_file_x=${work_dir}/${ref_ver}/${mt}/chrX.bed
GC_file_x=${work_dir}/${ref_ver}/${mt}/chrX.gc.bed
BED_file_y=${work_dir}/${ref_ver}/${mt}/chrY.bed
GC_file_y=${work_dir}/${ref_ver}/${mt}/chrY.gc.bed

### output files
out_prefix=${out_dir}/${ID}/${ID}
out_prefix_at=${out_prefix}_${region}
out_prefix_mt=${out_prefix}_chrM
out_prefix_mt_dup=${out_prefix}_chrM_dup
out_prefix_x=${out_prefix}_chrX
out_prefix_y=${out_prefix}_chrY

### make sure input file is bam or cram

if [[ "$input_fmt" != "bam" && "$input_fmt" != "cram" ]]
then
    echo 'Input file must be bam or cram file.'
    exit 1
fi

### coverage calculation
if [[ "$input_fmt" == "bam" ]]
then 
    $mosdepth -b $BED_file_at -n -x -F 3844 -Q 30 $out_prefix_at $BAM_file
    $mosdepth -b $BED_file_mt -n -x -F 3844 -Q 30 $out_prefix_mt $BAM_file
    $mosdepth -b $BED_file_mt -n -x -F 2820 -Q 30 $out_prefix_mt_dup $BAM_file
    $mosdepth -b $BED_file_x -n -x -F 3844 -Q 30 $out_prefix_x $BAM_file
    $mosdepth -b $BED_file_y -n -x -F 3844 -Q 30 $out_prefix_y $BAM_file
else
    $mosdepth -b $BED_file_at -f $ref_fasta -n -x -F 3844 -Q 30 $out_prefix_at $BAM_file
    $mosdepth -b $BED_file_mt -f $ref_fasta -n -x -F 3844 -Q 30 $out_prefix_mt $BAM_file
    $mosdepth -b $BED_file_mt -f $ref_fasta -n -x -F 2820 -Q 30 $out_prefix_mt_dup $BAM_file
    $mosdepth -b $BED_file_x -f $ref_fasta -n -x -F 3844 -Q 30 $out_prefix_x $BAM_file
    $mosdepth -b $BED_file_y -f $ref_fasta -n -x -F 3844 -Q 30 $out_prefix_y $BAM_file
fi

gunzip ${out_prefix_at}.regions.bed.gz
gunzip ${out_prefix_mt}.regions.bed.gz
gunzip ${out_prefix_mt_dup}.regions.bed.gz
gunzip ${out_prefix_x}.regions.bed.gz
gunzip ${out_prefix_y}.regions.bed.gz

paste <(grep -v '#' $GC_file_at) <(cut -f4 ${out_prefix_at}.regions.bed)|sed '1i chr\tstart\tend\tGC\tmap\tRC' > ${out_prefix_at}
paste <(grep -v '#' $GC_file_mt) <(cut -f4 ${out_prefix_mt}.regions.bed)|sed '1i chr\tstart\tend\tGC\tmap\tRC' > ${out_prefix_mt}
paste <(grep -v '#' $GC_file_mt) <(cut -f4 ${out_prefix_mt_dup}.regions.bed)|sed '1i chr\tstart\tend\tGC\tmap\tRC' > ${out_prefix_mt_dup}
paste <(grep -v '#' $GC_file_x) <(cut -f4 ${out_prefix_x}.regions.bed)|sed '1i chr\tstart\tend\tGC\tmap\tRC' > ${out_prefix_x}
paste <(grep -v '#' $GC_file_y) <(cut -f4 ${out_prefix_y}.regions.bed)|sed '1i chr\tstart\tend\tGC\tmap\tRC' > ${out_prefix_y}

rm ${out_prefix_at}.mosdepth.*
rm ${out_prefix_at}.regions.*
rm ${out_prefix_mt}.mosdepth.*
rm ${out_prefix_mt}.regions.*
rm ${out_prefix_mt_dup}.mosdepth.*
rm ${out_prefix_mt_dup}.regions.*
rm ${out_prefix_x}.mosdepth.*
rm ${out_prefix_x}.regions.*
rm ${out_prefix_y}.mosdepth.*
rm ${out_prefix_y}.regions.*

### mtDNA-CN estimation
Rscript ${work_dir}/mtDNA_CN.R ${out_prefix} ${region}

rm ${out_prefix_at}
rm ${out_prefix_mt}
rm ${out_prefix_mt_dup}
rm ${out_prefix_x}
rm ${out_prefix_y}

res_file=${out_prefix_at}.mitoCN.txt
if [[ -f "$res_file" ]]; then
    echo 'Complete successfully!'
fi
