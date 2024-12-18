# mitoCN

Estimate mitochondrial DNA copy number (mtDNA-CN) from whole-genome sequencing (WGS) data, with coverage bias adjustment, including GC bias and homology bias adjustment. It takes about 10 minutes of CPU time for a 30X genome. 

Please contact the author, Longfei Wang <wang.lo@wehi.edu.au>, if you would like to report any issues, feedback or feature requests.

## Citation
If you use mitoCN, please acknowledge by citing 
"Wang, L., Han, J., Fearnley, L.G. et al. Peripheral immune cell abundance differences link blood mitochondrial DNA copy number and Parkinson’s disease. npj Parkinsons Dis. 10, 219 (2024). https://doi.org/10.1038/s41531-024-00831-x"

## Installation
* Clone this repository
    ```
    git clone https://github.com/bahlolab/mitoCN
    ```
    
* Install [R](https://www.r-project.org/).

* Install [mosdepth](https://github.com/brentp/mosdepth). Add the mosdepth to PATH using 
    ```
    export PATH="/path/to/mosdepth:$PATH"
    ```


## Usage
* BAM input
    ```
    bash /path/to/mitoCN.sh -d <mitoCN_PATH> -f <BAM> -v <ref> -m <MT/chrM> -k <region> -o <out.prefix>
    ```
* CRAM input
    ```
    bash /path/to/mitoCN.sh -d <mitoCN_PATH> -f <CRAM> -a <FASTA> -v <ref> -m <MT/chrM> -k <region> -o <out.prefix>
    ```


**Params**  
* `-d` - the location of mitoCN. For example, if you download mitoCN to your home directory, the parameter should be `-d /Home/mitoCN`.
* `-f` - BAM/CRAM file name, must have index file. Single file only.
* `-v` - hg19, hg38. Depends on which reference you used for alignment.
* `-m` - “MT”: 1,2,…,24,MT; “chrM”: chr1,chr2,…,chrY,chrM. You can check this from the header using `samtools view -H $BAM > header`.
* `-k` - chr1, chr20, k500, by default k500 (recommend!). `-k k500` means select 500 sets of read bins with the same GC content distribution of mtDNA from each autosome. 
* `-o` - output file directory.
* `-a` - reference (.fasta) file name, if input file is CRAM format. Must include the index (.fai) file in the same directory.


## Example
Simulated reads were generated using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm).
Fold coverage = 10X, gender = male (1 copy of chrX & chrY), true mtDNA copy number = 4. The example data can be downloaded [here](https://zenodo.org/record/7964357). The example data should be put into the example folder and named as `data.tar.gz`.
The script of generating simulation data can be found [here](https://github.com/bahlolab/mitoCN/blob/main/example/reads_sim.sh).

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

## Run time
Use a GTEx sample as an example. ID = GTEX-1117F-0003-SM-6WBT7, CRAM file size = 46.3GB, fold coverage = 30X. The run time with difference -k options are as follows.

option | real time | user time | system time
--- | --- | --- | ---
k500 | 10m23s | 5m12s | 5m7s
chr1 | 1m50s | 56s | 50s
chr20 | 59s | 29s | 27s

## Docker Usage

### Running
A pre-built docker image is provided for convenience.
Ensure all your input data (crams, fastas etc are in a single `data` directory), then:
```bash
docker run --volume /path/to/data:/data ghcr.io/wehi-researchcomputing/mitocn -f /data/sample.cram -a /data/reference.fa -v hg38 -m chrM -k k500 -o /tmp/out
```

### Building 
Note, this isn't needed for normal usage, only for development or to make modifications to the image.
```bash
git clone git@github.com:bahlolab/mitoCN.git
docker build mitoCN --tag ghcr.io/wehi-researchcomputing/mitocn
```

## WDL Usage

`wf.wdl` provides a way to run mitoCN as part of a WDL workflow, which can allow you to use it in cloud environments such as Terra.
The WDL task pulls the premade docker image, meaning that simply uploading the `wf.wdl` is all you need to enable this workflow.

## Funding
This work is supported by the Michael J. Fox Foundation and the Shake It Up Australia Foundation through Spring 2022 RFA: Parkinson’s Pathway Molecular Data Analysis Program. 
