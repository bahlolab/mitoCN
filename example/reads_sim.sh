# ART-MountRainier-2016-06-05
# art_illumina -ss HSXn -sam -i <seq_ref_file> -l 150 -f <fold_coverage> -o <outfile_prefix>

# extract chr20, chrX, chrY, chrM from reference fa file
samtools faidx Homo_sapiens_assembly38.fasta chrM > ./reads_sim/human38_chrM.fa
samtools faidx Homo_sapiens_assembly38.fasta chrX > ./reads_sim/human38_chrX.fa
samtools faidx Homo_sapiens_assembly38.fasta chrY > ./reads_sim/human38_chrY.fa
samtools faidx Homo_sapiens_assembly38.fasta chr20 > ./reads_sim/human38_chr20.fa
samtools faidx Homo_sapiens_assembly38.fasta chr20 chrX chrY chrM > ./reads_sim/human38_select.fa

# simulate reads: male, mtDNA-CN = 3, average coverage = 6X
cd reads_sim
/home/users/allstaff/wang.lo/lab_bahlo_home/software/art_bin_MountRainier/art_illumina -ss HSXn -sam -i human38_chrM.fa -p -l 150 -f 9 -m 200 -s 10 -o chrM
/home/users/allstaff/wang.lo/lab_bahlo_home/software/art_bin_MountRainier/art_illumina -ss HSXn -sam -i human38_chrX.fa -p -l 150 -f 3 -m 200 -s 10 -o chrX
/home/users/allstaff/wang.lo/lab_bahlo_home/software/art_bin_MountRainier/art_illumina -ss HSXn -sam -i human38_chrY.fa -p -l 150 -f 3 -m 200 -s 10 -o chrY
/home/users/allstaff/wang.lo/lab_bahlo_home/software/art_bin_MountRainier/art_illumina -ss HSXn -sam -i human38_chr20.fa -p -l 150 -f 6 -m 200 -s 10 -o chr20

# sort and index sam -> bam
samtools sort chrM.sam > chrM.bam
samtools sort chrX.sam > chrX.bam
samtools sort chrY.sam > chrY.bam
samtools sort chr20.sam > chr20.bam

samtools index chrM.bam
samtools index chrX.bam
samtools index chrY.bam
samtools index chr20.bam

# merge bam files
samtools merge sample.bam chr20.bam chrX.bam chrY.bam chrM.bam
samtools index sample.bam

# convert bam to cram
samtools view -T human38_select.fa -C -o sample.cram sample.bam
samtools index sample.cram





