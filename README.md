# mitoCN

bash mitoCN.sh -d <mitoCN.dir> -f <BAM/CRAM> -v <ref> -m <MT/chrM> -k <region> -o <out.prefix>

-d: software directory
-f: BAM/CRAM file name, must have index file
-v: hg19, hg38
-m: “MT”: 1,2,…,24,MT; “chrM”: chr1,chr2,…,chrY,chrM
-k: auto, chr1, chr21, k500, by default k500
-o: output file location and prefix name
-a: fasta file, if input file is CRAM format
