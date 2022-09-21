#!/usr/bin/env nextflow
/*================================================================

								---- MODULE PIPELINE ---------

/*================================================================
The Aguilar Lab presents...

- A pipeline to extract and create miRNA and 3'UTR consensus sequences for analysis
   with targetscan and miRmap.

==================================================================
Version: 0.1
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

process COMPARE_TARGETS_TOOLS {
	tag "$CHR"
	echo true

	publishDir "${intermediates_dir}/compare-tools/",mode:"symlink"

	input:
	tuple val(CHR), file(TSOUT), file(MIRMAP)
	each BED
  each Rscript

	output:
	file "*.png"
	tuple val(CHR), path("*.tsv"), emit: TSV

	shell:
	"""
	echo " \$USER"
	echo " \$PATH"
  Rscript --vanilla ${Rscript} ${TSOUT} ${MIRMAP} ${BED} ${CHR}${params.output_name}

	"""
	stub:
	"""
	      touch  targets.${params.output_name}.tsv
				touch  targets.${params.output_name}.png
	"""
}
