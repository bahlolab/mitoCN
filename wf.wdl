version 1.0

struct IndexedFile {
    File main_file
    File index
}

task CombinedAnalysis {
    input {
        # See below for parameter descriptions
        IndexedFile cram
        String sample_id
        IndexedFile reference
        File mosdepth_bed
        File mosdepth_mt_bed

        # Sensible default resources
        String memory = "25G"
        Int cpus = 1
        String disk_space = "local-disk ${disk_space_int} HDD"
    } 
    parameter_meta {
        cram: {
            help: "Alignment CRAM to process."
        }
        sample_id: {
            help: "String uniquely identifying the sample. The output file will be named using this identifier so it must be unique."
        }
        reference: {
            help: "Reference genome in the fasta format, including an .fai index."
        }
        mosdepth_bed: {
            help: "BED file to be used to calculate the de-duplicated read depth statistics."
        }
        mosdepth_mt_bed: {
            help: "BED file to be used to target the mitochondria to calculate duplicated mitochondrial read depth statistics."
        }
        memory: {
            help: "Amount of RAM to request for this job, as a string, as described in https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#memory."
        }
        cpus: {
            help: "Optional integer indicating the number of CPUs to allocate."
        }
        disk_space: {
            help: "An optional disk space specification. If provided, should be as described in https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#disks."
        }   
    }
    # The default disk space is 120% of the input files, rounded up to the nearest integer
    Int disk_space_int = ceil(1.2 * size([
        reference.main_file,
        reference.index,
        mosdepth_bed,
        cram.main_file,
        cram.index
    ], "G"))
    runtime {
       docker: "quay.io/biocontainers/mosdepth:0.3.3--hd299d5a_3"
       disks: disk_space
       cpu: cpus
       memory: memory
    }
    command {
        set -e

        mosdepth \
            --by ${mosdepth_bed} \
            --fasta ${reference.main_file} \
            --no-per-base \
            --fast-mode \
            --flag 3844 \
            --quantize 30 \
            mosdepth-results \
            ${cram.main_file} &

        mosdepth \
            --by ${mosdepth_mt_bed} \
            --fasta ${reference.main_file} \
            --no-per-base \
            --fast-mode \
            --flag 2820 \
            --quantize 30 \
            mosdepth-results-dup \
            ${cram.main_file} &

        wait
    }
    output {
        Array[File] regions = glob("mosdepth-results.regions*")
        Array[File] dup_regions = glob("mosdepth-results-dup.regions*")
    }
}

workflow wf {
    # The task inputs implicitly become workflow inputs
    call CombinedAnalysis

    # Explicitly defining workflow outputs isn't necessary in WDL, but it is neccessary to use Terra's 
    # data table mapping feature
    output {
        Array[File] regions = CombinedAnalysis.regions
        Array[File] dup_regions = CombinedAnalysis.dup_regions
    }
}
