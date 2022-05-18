#!/usr/bin/env nextflow
/*================================================================

								---- MODULE PIPELINE ---------

/*================================================================
The Aguilar Lab presents...

- A pipeline to classify SNPs in microRNA regions and provide an overview of
diseases associated with microRNAs that present SNPs

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


process CAT_TARGETSCAN {
	tag "$TSOUT"

	publishDir "${results_dir}/cat-targets/",mode:"copy"

	input:
	file TSOUT

	output:
	file "${params.output_name}"

	shell:
  """
	cat *.tsout \
	| grep -P -v "a_Gene_ID\tmiRNA_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\tGroup_num\tSite_type\tmiRNA in this species\tGroup_type\tSpecies_in_this_group\tSpecies_in_this_group_with_this_site_type\tORF_overlap" > all.tmp
	echo "GeneID\tmiRNA_ID\tspecies_ID\tMSA_start\tMSA_end\tUTR_start\tUTR_end\tGroup_num\tSite_type\tmiRNA in this species\tGroup_type\tSpecies_in_this_group\tSpecies_in_this_group_with_this_site_type\tORF_overlap" > header.txt
	cat header.txt all.tmp | cut -f 1,2,6,7,9 > ${params.output_name}

	"""
	stub:
	"""
	     touch ${params.output_name}
	"""
}
