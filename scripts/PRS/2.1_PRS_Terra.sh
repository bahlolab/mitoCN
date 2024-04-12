# Install [PlINK 2.0](https://www.cog-genomics.org/plink/2.0/)
# 32 CPUs, 120G REM, 500GB disk
cd /home/jupyter/software
wget https://s3.amazonaws.com/plink2-assets/alpha4/plink2_linux_x86_64_20230829.zip
unzip plink2_linux_x86_64_20230829.zip
rm plink2_linux_x86_64_20230829.zip

work_dir="/home/jupyter/PRS"

mkdir ${work_dir}
cd ${work_dir}

plink="/home/jupyter/software/plink2"

# copy plink files
# gsutil -u terra-3f9b7ae9 cp gs://amp-pd-genomics/releases/2022_v3release_1115/wgs-WB-DWGS/plink/pfiles/all_chrs_merged.p* ./data/
gsutil -u terra-3f9b7ae9 cp gs://amp-pd-genomics/releases/2021_v2-5release_0510/wgs/plink/bfiles/all_chrs_merged* .

plink_file=${work_dir}/all_chrs_merged

### GRS file
### rsID  ALT-ALLELE  BETA
gsutil cp 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/adj_mtCN_GRS.txt' adj_mtCN_GRS.txt
gsutil cp 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/raw_mtCN_GRS.txt' raw_mtCN_GRS.txt
# head META5_GRS_RSid.txt
head adj_mtCN_GRS.txt
head raw_mtCN_GRS.txt

# GRS_path=${work_dir}/META5_GRS_RSid.txt
adj_mtCN_GRS_path=${work_dir}/adj_mtCN_GRS.txt
raw_mtCN_GRS_path=${work_dir}/raw_mtCN_GRS.txt

# out_path=${work_dir}/PD_GRS_loci
adj_mtCN_out_path=${work_dir}/adj_mtCN_GRS_loci
raw_mtCN_out_path=${work_dir}/raw_mtCN_GRS_loci

# extract the variants from plink files
# $plink --bfile ${plink_file} --extract ${GRS_path} --make-bed --out ${out_path}
$plink --bfile ${plink_file} --extract ${adj_mtCN_GRS_path} --make-bed --out ${adj_mtCN_out_path}
$plink --bfile ${plink_file} --extract ${raw_mtCN_GRS_path} --make-bed --out ${raw_mtCN_out_path}

# !check how many variants extracted
# cat ${out_path}.log
cat ${adj_mtCN_out_path}.log #84
cat ${raw_mtCN_out_path}.log #123

gsutil cp adj_mtCN_GRS_loci.bim 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/'
gsutil cp raw_mtCN_GRS_loci.bim 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/'

gsutil cp 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/missing_mtcn.txt' missing_adj_mtcn.txt
$plink --bfile ${plink_file} --extract range ${work_dir}/missing_adj_mtcn.txt --make-bed --out ${work_dir}/temp_adj
gsutil cp temp_adj.bim 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/'
$plink --bfile ${plink_file} --chr 10 --from-kb 58393 --to-kb 58394 --make-bed --out ${work_dir}/temp

gsutil cp 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/missing_raw_mtcn.txt' missing_raw_mtcn.txt
$plink --bfile ${plink_file} --extract range ${work_dir}/missing_raw_mtcn.txt --make-bed --out ${work_dir}/temp_raw
gsutil cp temp_raw.bim 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/'
$plink --bfile ${plink_file} --chr X --from-kb 46735 --to-kb 46736 --make-bed --out ${work_dir}/temp_raw

# replace the SNPs for PRS
gsutil cp 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/adj_mtCN_GRS_replace5snps.txt' adj_mtCN_GRS_replace5snps.txt
gsutil cp 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/raw_mtCN_GRS_replace11snps.txt' raw_mtCN_GRS_replace11snps.txt

adj_mtCN_GRS_replace_path=${work_dir}/adj_mtCN_GRS_replace5snps.txt
raw_mtCN_GRS_replace_path=${work_dir}/raw_mtCN_GRS_replace11snps.txt

# extract variants again
$plink --bfile ${plink_file} --extract ${adj_mtCN_GRS_replace_path} --make-bed --out ${adj_mtCN_out_path}
$plink --bfile ${plink_file} --extract ${raw_mtCN_GRS_replace_path} --make-bed --out ${raw_mtCN_out_path}
# !check how many variants extracted
cat ${adj_mtCN_out_path}.log
cat ${raw_mtCN_out_path}.log

# calculate PRS
$plink --bfile ${adj_mtCN_out_path} --score ${adj_mtCN_GRS_replace_path} list-variants --out ${work_dir}/adj_mtcn_PRS
gsutil cp adj_mtcn_PRS.sscore 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/'

$plink --bfile ${raw_mtCN_out_path} --score ${raw_mtCN_GRS_replace_path} list-variants --out ${work_dir}/raw_mtCN_PRS
gsutil cp raw_mtCN_PRS.sscore 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PRS/'