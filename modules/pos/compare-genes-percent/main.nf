#!/usr/bin/env nextflow
/*================================================================

								---- MODULE PIPELINE ---------

/*================================================================
The Aguilar Lab presents...

- A pipeline to extract and create miRNA and 3'UTR consensus sequences for analysis
   with targetscan and miRmap.

==================================================================
Version: 0.2
Project repository:
==================================================================
Authors:

- Bioinformatics Design
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)



- Bioinformatics Development
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)


- Nextflow Port
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)

///////////////////////////////////////////////////////////////

  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf-compare-miRNome"

/*This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*================================================================/*

/* MODULE START */

/* PRE1_CONVERT_GFF_TO_BED */

process COMPARE_GENES_PERCENT {
	tag "$REF, $ALT"

	publishDir "${results_dir}/compare-genes/",mode:"copy"

	input:
	file (REF)
	file (ALT)
  each Rscript

	output:
	file "*.png"
	path "*percent_filtered.tsv", emit: FILTERED_GENES
	path "*percent.tsv", emit: UNFILTERED_GENES
	file "*"

	shell:
	"""
  Rscript --vanilla  /usr/local/bin/compare_genes_percent.r ${REF} ${ALT} compare_genes_percent

	"""
	stub:
	"""
				touch  ${CHR}.changes.png
	"""
}
