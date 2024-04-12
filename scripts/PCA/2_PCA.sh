mkdir PCA
cd PCA

# Install [PlINK 2.0](https://www.cog-genomics.org/plink/2.0/)
# 32 CPUs, 120G REM, 500GB disk
uname -a
wget https://s3.amazonaws.com/plink2-assets/alpha4/plink2_linux_x86_64_20230829.zip
unzip plink2_linux_x86_64_20230829.zip
rm plink2_linux_x86_64_20230829.zip

# Download Data

## plink pfiles
mkdir data
gsutil -u terra-3f9b7ae9 cp gs://amp-pd-genomics/releases/2022_v3release_1115/wgs-WB-DWGS/plink/pfiles/all_chrs_merged.p* ./data/

## sample list
mkdir sample_list
gsutil cp 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/*_ID.txt' ./sample_list/

# QC
plink2='/home/jupyter/PCA/plink2'
pfile="/home/jupyter/PCA/data/all_chrs_merged"
afterQC="/home/jupyter/PCA/data/afterQC"

$plink2 --pfile $pfile --chr-set 29 --autosome --max-alleles 2 --geno 0.05 --maf 0.1 --make-bed --out $afterQC
# 10,418 samples
# 140714187 out of 158715874 variants
# 565619 variants removed due to missing genotype data
# 135504324 variants removed due to allele frequency threshold
# 4644244 variants remaining

mkdir results
rm ${pfile}*

# Run PCA analysis

## BioFIND
cohort='BioFIND'
### extract BioFIND
$plink2 --bfile $afterQC --keep ./sample_list/${cohort}_ID.txt --make-bed --out ./data/${cohort}
### QC within cohort
$plink2 --bfile ./data/${cohort} --chr-set 29 --autosome --max-alleles 2 --mind 0.1 --geno 0.05 --maf 0.1 --make-bed --out ./data/${cohort}_QC
### LD prunning
$plink2 -bfile ./data/${cohort}_QC --indep-pairwise 100kb 0.8 --seed 100 --out ./data/${cohort} --set-missing-var-ids @:#[b38]
### relatedness
$plink2 -bfile ./data/${cohort}_QC --king-cutoff 0.0884 --out ./data/${cohort}
### LD pruning + relatedness removal
$plink2 --bfile ./data/${cohort}_QC --extract ./data/${cohort}.prune.in --remove ./data/${cohort}.king.cutoff.out.id --make-bed --out ./data/${cohort}_pruned --set-missing-var-ids @:#[b38]
### run PCA
$plink2 --bfile ./data/${cohort}_pruned --chr-set 29 --pca --out ./results/plinkPCA_${cohort}

gsutil cp ./results/plinkPCA_${cohort}* 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'


## SURE
cohort='SURE'
### extract BioFIND
$plink2 --bfile $afterQC --keep ./sample_list/${cohort}_ID.txt --make-bed --out ./data/${cohort}
### QC within cohort
$plink2 --bfile ./data/${cohort} --chr-set 29 --autosome --max-alleles 2 --mind 0.1 --geno 0.05 --maf 0.1 --make-bed --out ./data/${cohort}_QC
### LD prunning
$plink2 -bfile ./data/${cohort}_QC --indep-pairwise 100kb 0.8 --seed 100 --out ./data/${cohort} --set-missing-var-ids @:#[b38]
### relatedness
$plink2 -bfile ./data/${cohort}_QC --king-cutoff 0.0884 --out ./data/${cohort}
### LD pruning + relatedness removal
$plink2 --bfile ./data/${cohort}_QC --extract ./data/${cohort}.prune.in --remove ./data/${cohort}.king.cutoff.out.id --make-bed --out ./data/${cohort}_pruned --set-missing-var-ids @:#[b38]
### run PCA
$plink2 --bfile ./data/${cohort}_pruned --chr-set 29 --pca --out ./results/plinkPCA_${cohort}

gsutil cp ./results/plinkPCA_${cohort}* 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'

## STEADY
cohort='STEADY'
### extract BioFIND
$plink2 --bfile $afterQC --keep ./sample_list/${cohort}_ID.txt --make-bed --out ./data/${cohort}
### QC within cohort
$plink2 --bfile ./data/${cohort} --chr-set 29 --autosome --max-alleles 2 --mind 0.1 --geno 0.05 --maf 0.1 --make-bed --out ./data/${cohort}_QC
### LD prunning
$plink2 -bfile ./data/${cohort}_QC --indep-pairwise 100kb 0.8 --seed 100 --out ./data/${cohort} --set-missing-var-ids @:#[b38]
### relatedness
$plink2 -bfile ./data/${cohort}_QC --king-cutoff 0.0884 --out ./data/${cohort}
### LD pruning + relatedness removal
$plink2 --bfile ./data/${cohort}_QC --extract ./data/${cohort}.prune.in --remove ./data/${cohort}.king.cutoff.out.id --make-bed --out ./data/${cohort}_pruned --set-missing-var-ids @:#[b38]
### run PCA
$plink2 --bfile ./data/${cohort}_pruned --chr-set 29 --pca --out ./results/plinkPCA_${cohort}

gsutil cp ./results/plinkPCA_${cohort}* 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'

## PPMI
cohort='PPMI'
### extract BioFIND
$plink2 --bfile $afterQC --keep ./sample_list/${cohort}_ID.txt --make-bed --out ./data/${cohort}
### QC within cohort
$plink2 --bfile ./data/${cohort} --chr-set 29 --autosome --max-alleles 2 --mind 0.1 --geno 0.05 --maf 0.1 --make-bed --out ./data/${cohort}_QC
### LD prunning
$plink2 -bfile ./data/${cohort}_QC --indep-pairwise 100kb 0.8 --seed 100 --out ./data/${cohort} --set-missing-var-ids @:#[b38]
### relatedness
$plink2 -bfile ./data/${cohort}_QC --king-cutoff 0.0884 --out ./data/${cohort}
### LD pruning + relatedness removal
$plink2 --bfile ./data/${cohort}_QC --extract ./data/${cohort}.prune.in --remove ./data/${cohort}.king.cutoff.out.id --make-bed --out ./data/${cohort}_pruned --set-missing-var-ids @:#[b38]
### run PCA
$plink2 --bfile ./data/${cohort}_pruned --chr-set 29 --pca --out ./results/plinkPCA_${cohort}

gsutil cp ./data/${cohort}.king.cutoff.out.id 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'
gsutil cp ./results/plinkPCA_${cohort}* 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'

## PDBP
cohort='PDBP'
### extract BioFIND
$plink2 --bfile $afterQC --keep ./sample_list/${cohort}_ID.txt --make-bed --out ./data/${cohort}
### QC within cohort
$plink2 --bfile ./data/${cohort} --chr-set 29 --autosome --max-alleles 2 --mind 0.1 --geno 0.05 --maf 0.1 --make-bed --out ./data/${cohort}_QC
### LD prunning
$plink2 -bfile ./data/${cohort}_QC --indep-pairwise 100kb 0.8 --seed 100 --out ./data/${cohort} --set-missing-var-ids @:#[b38]
### relatedness
$plink2 -bfile ./data/${cohort}_QC --king-cutoff 0.0884 --out ./data/${cohort}
### LD pruning + relatedness removal
$plink2 --bfile ./data/${cohort}_QC --extract ./data/${cohort}.prune.in --remove ./data/${cohort}.king.cutoff.out.id --make-bed --out ./data/${cohort}_pruned --set-missing-var-ids @:#[b38]
### run PCA
$plink2 --bfile ./data/${cohort}_pruned --chr-set 29 --pca --out ./results/plinkPCA_${cohort}

gsutil cp ./data/${cohort}.king.cutoff.out.id 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'
gsutil cp ./results/plinkPCA_${cohort}* 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'

## HBS
cohort='HBS'
### extract BioFIND
$plink2 --bfile $afterQC --keep ./sample_list/${cohort}_ID.txt --make-bed --out ./data/${cohort}
### QC within cohort
$plink2 --bfile ./data/${cohort} --chr-set 29 --autosome --max-alleles 2 --mind 0.1 --geno 0.05 --maf 0.1 --make-bed --out ./data/${cohort}_QC
### LD prunning
$plink2 -bfile ./data/${cohort}_QC --indep-pairwise 100kb 0.8 --seed 100 --out ./data/${cohort} --set-missing-var-ids @:#[b38]
### relatedness
$plink2 -bfile ./data/${cohort}_QC --king-cutoff 0.0884 --out ./data/${cohort}
### LD pruning + relatedness removal
$plink2 --bfile ./data/${cohort}_QC --extract ./data/${cohort}.prune.in --remove ./data/${cohort}.king.cutoff.out.id --make-bed --out ./data/${cohort}_pruned --set-missing-var-ids @:#[b38]
### run PCA
$plink2 --bfile ./data/${cohort}_pruned --chr-set 29 --pca --out ./results/plinkPCA_${cohort}

gsutil cp ./data/${cohort}.king.cutoff.out.id 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'
gsutil cp ./results/plinkPCA_${cohort}* 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'

## LBD
cohort='LBD'
### extract BioFIND
$plink2 --bfile $afterQC --keep ./sample_list/${cohort}_ID.txt --make-bed --out ./data/${cohort}
### QC within cohort
$plink2 --bfile ./data/${cohort} --chr-set 29 --autosome --max-alleles 2 --mind 0.1 --geno 0.05 --maf 0.1 --make-bed --out ./data/${cohort}_QC
### LD prunning
$plink2 -bfile ./data/${cohort}_QC --indep-pairwise 100kb 0.8 --seed 100 --out ./data/${cohort} --set-missing-var-ids @:#[b38]
### relatedness
$plink2 -bfile ./data/${cohort}_QC --king-cutoff 0.0884 --out ./data/${cohort}
### LD pruning + relatedness removal
$plink2 --bfile ./data/${cohort}_QC --extract ./data/${cohort}.prune.in --remove ./data/${cohort}.king.cutoff.out.id --make-bed --out ./data/${cohort}_pruned --set-missing-var-ids @:#[b38]
### run PCA
$plink2 --bfile ./data/${cohort}_pruned --chr-set 29 --pca --out ./results/plinkPCA_${cohort}

gsutil cp ./data/${cohort}.king.cutoff.out.id 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'
gsutil cp ./results/plinkPCA_${cohort}* 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'

## all
### LD prunning
$plink2 -bfile $afterQC --indep-pairwise 100kb 0.8 --seed 100 --out ./data/all --set-missing-var-ids @:#[b38]
### relatedness
$plink2 -bfile $afterQC --king-cutoff 0.0884 --out ./data/all
### LD pruning + relatedness removal
$plink2 --bfile $afterQC --extract ./data/all.prune.in --remove ./data/all.king.cutoff.out.id --make-bed --out ./data/all_pruned --set-missing-var-ids @:#[b38]
### run PCA
$plink2 --bfile ./data/all_pruned --chr-set 29 --pca --out ./results/plinkPCA_all

gsutil cp ./data/all.king.cutoff.out.id 'gs://fc-90503228-8a34-41b1-9199-faa3873a8cd2/PCA/'