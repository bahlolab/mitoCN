# mitoCN

Estimate mitochondrial DNA copy number (mtDNA-CN) from whole-genome sequencing (WGS) data, with coverage bias adjustment, including GC bias and homology bias adjustment.

## Installation
* Clone this repositoty
    ```
    git clone https://github.com/bahlolab/mitoCN
    ```
    
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
* `-d` - the location of mitoCN. For example, if you download mitoCN to your home direcoty, mitoCN_PATH=/Home/mitoCN.
* `-f` - BAM/CRAM file name, must have index file.
* `-v` - hg19, hg38. Depends on which reference you used for alignment.
* `-m` - “MT”: 1,2,…,24,MT; “chrM”: chr1,chr2,…,chrY,chrM. 
* `-k` - chr1, chr20, k500, by default k500 (recommend!).
* `-o` - output file directory.
* `-a` - reference (.fasta) file name, if input file is CRAM format. Must include the index (.fai) file in the same directory.


## Example
Simulated reads were generated using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). 
Fold coverage = 10X, gender = male, true mtDNA copy number = 4. The script can be found [here](https://github.com/bahlolab/mitoCN/blob/main/example/reads_sim.sh).

## Output
[Output file](https://github.com/bahlolab/mitoCN/blob/main/example/sample/sample_chr20.mitoCN.txt) of [example.sh](https://github.com/bahlolab/mitoCN/blob/main/example/example.sh)

region | Np | mt | var.mt | mt.dup | var.mt.dup | chrX | chrY
--- | --- | --- | --- |--- |--- |--- |---
chr20 | 10 | 3.6 | 0.04 | 3.6 | 0.04 | 0.87 | 0.86

* `region` - which region option was used, chr1, chr20, k500.
* `Np` - average coverage 
* `mt` - mtDNA-CN estimate without duplicated reads on mtDNA.
* `var.mt` - variance of mtDNA-CN estimate without duplicates.
* `mt.dup` - mtDNA-CN estimate with duplicated reads on mtDNA.
* `var.mt.dup` - variance of mtDNA-CN estimate with duplicates.
* `chrX` - X chromosome copy number estimate.
* `chrY` - Y chromosome copy number estimate.
