#!/usr/bin/env nextflow
nextflow.enable.dsl=2  // Enable Nextflow DSL2 syntax

// =====================
// Parameters
// =====================

// Input file listing paths to *_mitocalc.tar.gz files (one per line)
// Active input path (TSV format, still treated as plain text list of tarballs)
params.input_paths = "/vast/scratch/users/$USER/MTDNA/20240423_UKBB_mitocalc_paths.tsv"

// Output directory for summary result
params.outputDir = "/vast/scratch/users/$USER/mtDNA"

// Convert outputDir to Nextflow file object
outputDir = file(params.outputDir)

process processFiles {

	module 'R'                        // Load the R module (for Rscript)

	errorStrategy 'retry'            // Retry this process on failure

	cpus = 1                         // Allocate 1 CPU
	memory = { 500 * task.attempt + ' MB' }   // Memory scales with retry attempts
	time = { 30 * task.attempt + ' min'}      // Timeout scales with retry attempts

	input:
	path (input_file)               // Input is a path to one or more tarballs

	output:
	path "*_summary.txt"            // Output: summary text file per batch (batch = 500 tarballs)

	shell:
	'''
	# Loop over each mitocalc tarball file
	for file in !{input_file}
	do
		# Extract tarball
		tar -xvf $file
		
		# Derive sample ID from filename
		sample=$(basename "$file" '_mitocalc.tar.gz')

		# Run R script to compute mtDNA-CN from *_at, *_mt, *_x, *_y files
		Rscript /vast/projects/bahlo_ukbiobank/app36610_HDA/mitochondria/mtDNA_ukbb/mtDNA_CN_edit.R $sample k500

		# Extract summary line from result
		info=$(sed -n '2p' *_k500.mitoCN.txt)

		# Write sample + result line to output file
		echo -e "$sample\t$info" >> results_summary.txt

		# Clean up intermediate files
		rm *_k500.mitoCN.txt
		rm *_at
		rm *_mt
		rm *_x
		rm *_y
	done

	# Clean up tarballs after processing
	rm *.tar.gz
	'''
}

process combineFiles {

	publishDir "$outputDir", mode: "copy"  // Copy final output to outputDir

	input:
	path '*_summary.txt'  // Input: all the partial summary files

	output:
	path "final.txt"      // Output: single concatenated file

	script:
	"""
	# Add header to final file
	echo -e "sample\tregion\tNp\tm_beta\tk\tmt\tvar.mt\tchrX\tchrY" > final.txt
	
	# Append all individual summaries
	cat *_summary.txt >> final.txt

	# Clean up
	rm *_summary.txt
	"""
}

workflow {
	// Read paths from the input file and split them into chunks of 500 lines
	input_files = Channel
				.fromPath(params.input_paths)   // Read the file path
				.splitText()                    // Split lines (tarball paths)
				.collate( 500 )                 // Group into batches of 500 for efficient processing

	// Run processing and collect outputs
	processFiles(input_files)
	| collect         // Collect all outputs from batches
	| combineFiles    // Combine them into a final summary
}

