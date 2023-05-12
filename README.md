# mitoCN

Estimate mitochondrial DNA copy number (mtDNA-CN) from whole-genome sequencing (WGS) data, with coverage bias adjustment, including GC bias and homology bias adjustment.

## Installation
* Clone this repositoty
* Install [R](https://www.r-project.org/).
* Install [mosdepth](https://github.com/brentp/mosdepth). Add the location of `mosdepth` to mitoCN.sh using 
    ```
    mosdepth=/path/to/mosdepth/mosdepth_version
    ```

## Usage
BAM input

    ```
    bash /path/to/mitoCN.sh -d <mitoCN_PATH> -f <BAM> -v <ref> -m <MT/chrM> -k <region> -o <out.prefix>
    ```
CRAM input

    ```
    bash /path/to/mitoCN.sh -d <mitoCN_PATH> -f <CRAM> -a <FASTA> -v <ref> -m <MT/chrM> -k <region> -o <out.prefix>
    ```

**Params**  
* `-d` - path to mitoCH.sh.
* `-f` - BAM/CRAM file name, must have index file.
* `-v` - hg19, hg38.
* `-m` - “MT”: 1,2,…,24,MT; “chrM”: chr1,chr2,…,chrY,chrM.
* `-k` - chr1, chr20, k500, by default k500.
* `-o` - output file directory.
* `-a` - fasta file, if input file is CRAM format.

## Output


region | #Np | #mt | #var.mt | #mt.dup | #var.mt.dup | #chrX | #chrY
--- | --- | --- | --- |--- |--- |--- |---
k500 | 30.3 | 215.54 | 0.72 | 316.96 | 1.06 | 2.07 | 0
