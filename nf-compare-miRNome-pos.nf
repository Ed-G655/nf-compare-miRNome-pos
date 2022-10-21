#!/usr/bin/env nextflow

/*================================================================
The Aguilar Lab presents...

- A pipeline to compare microRNA targets from targetScan and miRmap
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

Pre-processing:

Core-processing:


Pos-processing
Pos1-COMPARE_TARGETS

Analysis:

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
	The miRNome compare pipeline
  v${version}
  ==========================================

	Usage:

	nextflow run ${pipeline_name}.nf --ref_dir <path to input 1> --alt_dir <path to input 2>
  --bed <path to input 3>    [--output_dir path to results ]

		--ref_dir	<- Directory whith the mirmap files;

	  --alt_dir	<- Directory whith the targetScan files;

    --bed <- bed file;

	  --output_dir     <- directory where results, intermediate and log files will be stored;
	      default: same dir where --query_fasta resides

	  -resume	   <- Use cached results if the executed project has been run before;
	      default: not activated
	      This native NF option checks if anything has changed from a previous pipeline execution.
	      Then, it resumes the run from the last successful stage.
	      i.e. If for some reason your previous run got interrupted,
	      running the -resume option will take it from the last successful pipeline stage
	      instead of starting over
	      Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf-miRNome-compare-pos"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.ref_dir = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.alt_dir = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.bed = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "${pipeline_name} v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.01.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
*/

/* Check if the input directory is provided
    if it was not provided, it keeps the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.ref_dir | !params.alt_dir | !params.bed ) {
  log.error " Please provide the --ref_dir AND --alt_dir AND --bed \n\n" +
  " For more information, execute: nextflow run ${pipeline_name} --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.ref_dir).getParent() //!! maybe creates bug, should check

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable pipeline_name defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
resulalt_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
Useful functions definition
*/

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The nf-miRNome-compare pos analysis pipeline
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['Input REF dir']			= params.ref_dir
pipelinesummary['Input ALT dir']			= params.alt_dir
pipelinesummary['Input bed']			= params.bed
pipelinesummary['Results Dir']		= resulalt_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/* enable DSL2*/
nextflow.enable.dsl=2

/*
	READ GENERAL INPUTS
*/

// Load mirmap REF data
Channel
	.fromPath( "${params.ref_dir}*.mirmapout" )
	.set{ mirmap_ref_input }

// Load mirmap alt data
	Channel
		.fromPath( "${params.alt_dir}*.alt.mirmapout" )
		.set{ mirmap_alt_input }

//load targetScan REF data
Channel
	.fromPath( "${params.ref_dir}*.tsout" )
	.set{ ts_ref_input }


	//load targetScan REF data
	Channel
		.fromPath( "${params.alt_dir}*.alt.tsout" )
		.set{ ts_alt_input }

// load BED file
Channel
	.fromPath( "${params.bed}" )
	.set{ bed_input }


/*
 Load R fileS
*/

/* R_script_5*/
Channel
	.fromPath( "./modules/pos/compare-tools/compare_tools.r" )
	.set{ R_script_5}

	/* R_script_6*/
	Channel
		.fromPath( "./modules/pos/compare-targets/compare_targets.r" )
		.set{ R_script_6}

	/* R_script_7*/
	Channel
					.fromPath( "./modules/pos/compare-mirnome/compare_mirnome.R" )
					.set{ R_script_7}

/* R_script_7*/
Channel
			.fromPath( "./modules/pos/eulerr-tools/euler.R" )
			.set{ R_script_8}

/* R_script_7*/
Channel
			 .fromPath( "./modules/pos/resume-changes/venn_data.R" )
			 .set{ R_script_9}
/* R_script_10*/
Channel
			 .fromPath( "./modules/pos/compare-targets_percent/compare_targets_percent.r" )
			 .set{ R_script_10}

/* Python script */
Channel
			.fromPath("./modules/pos/Venn_plot/compare_mirnome.py")
			.set{Python_script}
/* R_script_11*/
Channel
			 .fromPath( "./modules/pos/compare-genes-percent/compare_genes_percent.r" )
			 .set{ R_script_11}

/* R_script_12*/
Channel
			 .fromPath( "./modules/pos/compare-mirnome_genes/compare_mirnome_genes.R" )
			 .set{ R_script_12}

/* R_script_13*/
Channel
			 .fromPath( "./modules/pos/resume-changes-genes/venn_data_genes.R" )
			 .set{ R_script_13}

/* R_script_14*/
Channel
			 .fromPath( "./modules/pos/plot-filtered-genes/plot_filtered_genes.R" )
			 .set{ R_script_14}

/* R_script_15*/
Channel
			 .fromPath( "./modules/pos/plot-histogram-genes/plot_histogram_genes.R" )
			 .set{ R_script_15}

/* R_script_16*/
Channel
				.fromPath( "./modules/pos/resume-venn-genes/venn_data_genes.R" )
			 	.set{ R_script_16}



	/*	  Import modules */

								/*POS-processing */
include{CAT_TOOLS_TARGETS as CAT_REF_TARGETS} from './modules/pos/cat-tools-targets/main.nf' addParams(output_name: 'All_targets.ref.tsv')
include{CAT_TOOLS_TARGETS as CAT_ALT_TARGETS} from './modules/pos/cat-tools-targets/main.nf' addParams(output_name: 'All_targets.alt.tsv')


include{COMPARE_TARGETS_TOOLS as COMPARE_TOOLS_REF} from './modules/pos/compare-tools/main.nf' addParams(output_name: '.ref')
include{COMPARE_TARGETS_TOOLS as COMPARE_TOOLS_ALT} from './modules/pos/compare-tools/main.nf' addParams(output_name: '.alt')

include{COMPARE_TARGETS} from './modules/pos/compare-targets/main.nf'

include{COMPARE_TARGETS_PERCENT} from './modules/pos/compare-targets_percent/main.nf'

include{GREP_TARGETSID as GREP_TARGETSID_REF } from './modules/pos/grep-targets/main.nf' addParams(output_name: 'targetsID.ref.tsv')

include{GREP_TARGETSID as GREP_TARGETSID_ALT } from './modules/pos/grep-targets/main.nf' addParams(output_name: 'targetsID.alt.tsv')

include{GREP_TOOLS as GREP_TS_REF } from './modules/pos/grep-tools/main.nf' addParams(TOOLS: 'targetscan' , output_name: 'targets_TS.ref.tsv')
include{GREP_TOOLS as GREP_MP_REF } from './modules/pos/grep-tools/main.nf' addParams(TOOLS: 'mirmap' , output_name: 'targets_MP.ref.tsv')

include{GREP_TOOLS as GREP_TS_ALT } from './modules/pos/grep-tools/main.nf' addParams(TOOLS: 'targetscan' , output_name: 'targets_TS.alt.tsv')
include{GREP_TOOLS as GREP_MP_ALT } from './modules/pos/grep-tools/main.nf' addParams(TOOLS: 'mirmap' , output_name: 'targets_MP.alt.tsv')

include{COMPARE_MIRNOME} from './modules/pos/compare-mirnome/main.nf'

include{EULERR_MIRNOME} from './modules/pos/eulerr-tools/main.nf'

include{VENN_PLOT} from './modules/pos/Venn_plot/main.nf'

include{RESUME_CHANGES} from './modules/pos/resume-changes/main.nf'

include{RESUME_VENN_CHANGES} from './modules/pos/resume-venn-genes/main.nf'

include{COMPARE_GENES_PERCENT} from './modules/pos/compare-genes-percent/main.nf'

include{COMPARE_GENES_MIRNOME} from './modules/pos/compare-mirnome_genes/main.nf'

include{RESUME_GENE_CHANGES} from './modules/pos/resume-changes-genes/main.nf'

include{PLOT_FILTERED_GENES} from './modules/pos/plot-filtered-genes/main.nf'

include{CAT_PERCENT_GENES} from './modules/pos/cat-percents-genes/main.nf'

include{CAT_PERCENT_GENES as CAT_FILTERED_GENES} from './modules/pos/cat-percents-genes/main.nf'

include{PLOT_HISTOGRAM} from './modules/pos/plot-histogram-genes/main.nf'

/*  main pipeline logic */
workflow  {
/* pos-processing */
// Define function to get chrom
def get_chrom = { file -> file.baseName.replaceAll(/.alt/,"").replaceAll(/.filtered/,"")}

					// collect targets outputs
						// TARGETSCAN_REF = ts_ref_input.map{file -> tuple(file.baseName.replaceAll(/.filtered/,""), file) }
						TARGETSCAN_REF = ts_ref_input
						TARGETSCAN_ALT = ts_alt_input

						MIRMAP_REF = mirmap_ref_input
						MIRMAP_ALT = mirmap_alt_input


						// join channels
						REF_TARGETS_TOOLS = TARGETSCAN_REF.mix(MIRMAP_REF)
						ALT_TARGETS_TOOLS = TARGETSCAN_ALT.mix(MIRMAP_ALT)
						//
						// // Merge mirmap and targetscan data
						REF_TARGETS = COMPARE_TOOLS_REF(REF_TARGETS_TOOLS.collect(), bed_input, R_script_5)
					  // // Merge mirmap and targetscan data
					  ALT_TARGETS = COMPARE_TOOLS_ALT(ALT_TARGETS_TOOLS.collect(), bed_input, R_script_5)
						//
						//
						// // COMPARE_TARGETS: Compare REF and ALT targets
						COMPARE_TARGETS(REF_TARGETS.TSV, ALT_TARGETS.TSV, R_script_6)
						//
						COMPARE_TARGETS_PERCENT(REF_TARGETS.TSV, ALT_TARGETS.TSV, R_script_10)
						//CAT TARGETS OUTPUTS
						//CAT_TARGETS(COMPARE_TARGETS.out.CHANGES.collect())
						//
						// // CAT REF_TARGETS
						// CAT_REF_TARGETS(REF_TARGETS.TSV.collect())
						// // CAT ALT TARGETS
						// CAT_ALT_TARGETS(ALT_TARGETS.TSV.collect())
						//
						// // GREP REF TARGETSID
						// GREP_TARGETSID_REF(CAT_REF_TARGETS.out)
						// // GREP ALT TARGETSID
						// GREP_TARGETSID_ALT(CAT_ALT_TARGETS.out)

						// PLOT miRNome changes
						COMPARE_MIRNOME(REF_TARGETS.TSV, ALT_TARGETS.TSV, R_script_7)

						// PLOT TARGET TOOLS
			 			EULERR_MIRNOME(REF_TARGETS.TSV, ALT_TARGETS.TSV, R_script_8)
						//
						// // Sum changes data
						// RESUME_CHANGES(COMPARE_MIRNOME.out.VENN_DATA.collect(), R_script_9)
						// // Compare genes targets
						// COMPARE_GENES_PERCENT(REF_TARGETS.TSV, ALT_TARGETS.TSV, R_script_11)
						// // Compare miRNome
						// COMPARE_GENES_MIRNOME(REF_TARGETS.TSV, ALT_TARGETS.TSV, R_script_12) //
						// // Compare resume
						// RESUME_GENE_CHANGES(COMPARE_GENES_MIRNOME.out.VENN_DATA.collect(), R_script_13)
						//
						// RESUME_VENN_CHANGES(COMPARE_GENES_MIRNOME.out.VENN_DATA_GENES.collect(), R_script_16)
						// // // PLOT miRNome changes
						// // VENN_PLOT(CAT_REF_TARGETS.out, CAT_ALT_TARGETS.out, Python_script)
						// CAT_FILTERED_GENES(COMPARE_GENES_PERCENT.out.FILTERED_GENES.collect())
						// // PLot FILTERED_GENES
						// PLOT_FILTERED_GENES(CAT_FILTERED_GENES.out, R_script_14)
						// // CAT UNFILTERED GENES
						// CAT_PERCENT_GENES(COMPARE_GENES_PERCENT.out.UNFILTERED_GENES.collect())
						// // PLOT HISTOGRAM
					 	// PLOT_HISTOGRAM(CAT_PERCENT_GENES.out, R_script_15)
}
