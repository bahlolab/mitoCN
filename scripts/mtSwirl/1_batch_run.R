rm(list=ls())

#!!! study <- "LBD"
IDs <- list.dirs(paste0("/vast/scratch/users/wang.lo/mtSwirl/output/res/",study,"/"), full.names = F)[-1]
template <- readLines("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/mtSwirl/script/single_sample.sh")
template <- template[!(substring(template, 1, 1) %in% c("", "#"))]
work_dir <- "/vast/scratch/users/wang.lo/mtSwirl/run/"
if(!(dir.exists(work_dir))) dir.create(work_dir,recursive = T)
script_dir <- paste0(work_dir,"script/")
if(!(dir.exists(script_dir))) dir.create(script_dir,recursive = T)
run_sh <- paste0(work_dir,"run_",study,".sh")
cat(paste0("cd ",script_dir,"\n"), file=run_sh)

for(i in 1:length(IDs)){
  
  id <- IDs[i]
  print(id)
  
  sh <- paste0(script_dir,id,".sh")
  mem=4; cpu=1
  cat("#!/bin/bash\n\n", file=sh)
  cat(paste0("#SBATCH --time=24:00:00\n#SBATCH --mem=",mem,"G\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=",cpu,"\n#SBATCH --job-name=",id,"\n#SBATCH --error=",id,".err\n#SBATCH --output=",id,".out\n#SBATCH --mail-user=wang.lo@wehi.edu.au\n\n"), file=sh, append=TRUE)
  
  cat(paste0("ID=",id,"\n"), file=sh, append = TRUE)
  cat(paste0("study=",study,"\n"), file=sh, append = TRUE)
  cat(paste(template, collapse="\n"), file=sh, append=TRUE)
  
  cat(paste0("sbatch ",sh,"\n"), file=run_sh, append = TRUE)
  
}

system(paste0("bash ",work_dir,"run_",study,".sh"))
