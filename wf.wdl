version 1.0

struct IndexedFile {
    File main_file
    File index
}

task mitoCN {
    input {
        IndexedFile cram
        IndexedFile reference
        # See docs for the -v flag
        String reference_version
        # See docs for the -m flag
        String mt_name
        # See docs for the -k flag
        Int bins

        # Sensible default resources
        String memory = "25G"
        Int cpus = 1
        String disk_space = "local-disk ${disk_space_int} HDD"
    } 
    # The default disk space is 120% of the input files, rounded up to the nearest integer
    Int disk_space_int = ceil(1.2 * size([
        reference.main_file,
        reference.index,
        cram.main_file,
        cram.index
    ], "G"))
    runtime {
       docker: "ghcr.io/wehi-researchcomputing/mitocn"
       disks: disk_space
       cpu: cpus
       memory: memory
    }
    command {
        ./mitoCN.sh -d /app -f ~{cram.main_file} -a ~{reference.main_file} -v ~{reference_version} -m ~{mt_name} -k ~{bins} -o ./out
    }
    output {
        File result = glob("out/*.txt")[0]
    }
}

workflow wf {
    # The task inputs implicitly become workflow inputs
    call mitoCN

    # Explicitly defining workflow outputs isn't necessary in WDL, but it is neccessary to use Terra's 
    # data table mapping feature
    output {
        File result = mitoCN.result
    }
}