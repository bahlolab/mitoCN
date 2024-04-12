
# ID="BF-1001"
# study="BioFIND"

module load R
dir_res=/vast/scratch/users/wang.lo/mtSwirl/output/res/${study}
dir_script="/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/mtSwirl/script"

# process flagstat, idxstat, yield matrics
flagstat=${dir_res}/${ID}/${ID}.flagstat.txt
flagstat_tmp=${dir_res}/${ID}/${ID}.flagstat.tmp.txt
Rscript ${dir_script}/reformat.flagstat.R ${flagstat} ${flagstat_tmp}

idxstat=${dir_res}/${ID}/${ID}.stats.tsv
idxstat_tmp=${dir_res}/${ID}/${ID}.idxstat.txt
Rscript ${dir_script}/idxstat.summary.R ${idxstat} ${idxstat_tmp}

yield=${dir_res}/${ID}/${ID}.yield_metrics.txt
yield_tmp=${dir_res}/${ID}/${ID}.yield_metrics.tmp.txt
cat ${yield} | tail -n 4 > ${yield_tmp}

# merge three files
summary_stat=${dir_res}/${ID}/${ID}.stat.summary.txt
paste -d "\t" ${yield_tmp} ${flagstat_tmp} ${idxstat_tmp} > ${summary_stat}

# calculate nuc_mean_coverage and mtcn
stat=${dir_res}/${ID}/${ID}_mtanalysis_diagnostic_statistics.tsv
mtcn=${dir_res}/${ID}/${ID}.mtcn.txt
Rscript ${dir_script}/nuc_mean_coverage.R ${summary_stat} ${stat} ${mtcn}

# archive files
dir_arv=/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/mtSwirl/
cp ${mtcn} ${dir_arv}/mtcn/${study}/
cp ${dir_res}/${ID}/${ID}.self.ref.split.selfToRef.final.vcf ${dir_arv}/chrM_variants/${study}/
cp ${dir_res}/${ID}/${ID}.appended.liftedOver.tsv ${dir_arv}/chrM_base_coverage/${study}/
cp ${summary_stat} ${dir_arv}/summary_statistic/${study}/
