# mitoCN

Estimate mitochondrial DNA copy number (mtDNA-CN) from whole-genome sequencing (WGS) data, with coverage bias adjustment, including GC bias and homology bias adjustment.

## Installation
* Clone this repositoty
* Install [mosdepth](https://github.com/brentp/mosdepth). Set the path using 
    ```
    export PATH=$PATH:/path/to/mosdepth
    ```

## Usage
* BAM input.

    ```
    bash /path/to/mitoCN.sh -d <mitoCN_PATH> -f <BAM> -v <ref> -m <MT/chrM> -k <region> -o <out.prefix>
    ```
* CRAM input.

    ```
    bash /path/to/mitoCN.sh -d <mitoCN_PATH> -f <CRAM> -a <FASTA> -v <ref> -m <MT/chrM> -k <region> -o <out.prefix>
    ```

* **Params**  
  * `-d` - path to mitoCH.sh.
  * `-f` - BAM/CRAM file name, must have index file.
  * `-v` - hg19, hg38.
  * `-m` - “MT”: 1,2,…,24,MT; “chrM”: chr1,chr2,…,chrY,chrM.
  * `-k` - auto, chr1, chr21, k500, by default k500.
  * `-o` - output file location and prefix name.
  * `-a` - fasta file, if input file is CRAM format.

## In-house data application
Please email me the following information
* BAM/CRAM file directory
* fasta file location if CRAM files
* metadata including sex, age, phenotypes
* reference file: hg19 or hg38
