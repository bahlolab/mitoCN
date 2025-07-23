version 1.0

task retools {
  # Define inputs to the task
  input {
    Array[File]+ cramfiles     # Input CRAM files (mapped sequencing reads)
    Array[File]+ craifiles     # Input CRAI index files
    File reference             # Reference genome FASTA
    File cache                 # Pre-populated HTSlib cache (tar.gz)
    File payload               # Tarball containing BED and GC annotation files
  }

  # Shell commands to run in the container
  command <<<
    set -euxo pipefail    # Enable strict shell options for safer execution

    # Unpack BED/gc annotation files and HTSlib cache
    tar xzf ~{payload}
    tar xzf ~{cache}

    # Iterate over all CRAM files
    for CRAM_file in ~{sep=" " cramfiles}; do

        # Derive base sample name by removing known suffix
        fname=$(basename $CRAM_file "_23193_0_0.cram")

        # Define output prefixes for different chromosomes
        out_prefix_at="${fname}_at"   # autosomal (chr1-22)
        out_prefix_mt="${fname}_mt"   # chrM
        out_prefix_x="${fname}_x"     # chrX
        out_prefix_y="${fname}_y"     # chrY

        # Run mosdepth to calculate read depth for each region
        mosdepth -b k500.bed -f ~{reference} -n -x -F 3844 -Q 30 ${out_prefix_at} ${CRAM_file}
        mosdepth -b chrM.bed -f ~{reference} -n -x -F 3844 -Q 30 ${out_prefix_mt} ${CRAM_file}
        mosdepth -b chrX.bed -f ~{reference} -n -x -F 3844 -Q 30 ${out_prefix_x} ${CRAM_file}
        mosdepth -b chrY.bed -f ~{reference} -n -x -F 3844 -Q 30 ${out_prefix_y} ${CRAM_file}

        # Decompress the mosdepth output region files
        gunzip ${out_prefix_x}.regions.bed.gz
        gunzip ${out_prefix_y}.regions.bed.gz
        gunzip ${out_prefix_mt}.regions.bed.gz
        gunzip ${out_prefix_at}.regions.bed.gz

        # Combine GC content annotations with read count results for each region
        paste k500.gc.bed <(cut -f4 ${fname}_at.regions.bed) | sed '1i chr\tstart\tend\tGC\tmap\tRC' > ${fname}_at
        paste chrM.gc.bed <(cut -f4 ${fname}_mt.regions.bed) | sed '1i chr\tstart\tend\tGC\tmap\tRC' > ${fname}_mt
        paste chrX.gc.bed <(cut -f4 ${fname}_x.regions.bed) | sed '1i chr\tstart\tend\tGC\tmap\tRC' > ${fname}_x
        paste chrY.gc.bed <(cut -f4 ${fname}_y.regions.bed) | sed '1i chr\tstart\tend\tGC\tmap\tRC' > ${fname}_y

        # Bundle outputs for mitoCN
        tar czvf ${fname}_mitocalc.tar.gz \
          ${fname}_y \
          ${fname}_x \
          ${fname}_at \
          ${fname}_mt

    done
  >>>

  # Output the tarballs created for each sample
  output {
    Array[File] mitocalc = glob("*mitocalc.tar.gz")   # GC-normalized depth files
  }

  # Runtime configuration for execution environment
  runtime {
      docker: "dx://project-GbZQf58JK9yfGYG8vk7qp6g9:/dxfuse_test_0.1_tar.gz"  # Docker image path on DNAnexus
      dx_timeout: "8H"                             # Set max execution time to 8 hours
      dx_instance_type: "mem3_ssd1_v2_x2"          # Instance type for moderate memory/SSD
  }

  parameter_meta {
    cramfiles: {
      description: "mapped short read data",
      patterns: ["*.cram"],
      stream: true
    }
    craifiles: {
      description: "mapped short read data indices",
      patterns: ["*.crai"],
      stream: true
    }
    reference: {
      description: "reference genome",
      patterns: ["*.fa","*.fasta"]
    }
    cache: {
      description: "htslib cache (prepopulated, tar.gz)",
      patterns: ["*.tar.gz"]
    }
    payload: {
      description: "payload file containing mitoCN BED files",
      patterns: ["*.tar.gz"]
    }
  }
}