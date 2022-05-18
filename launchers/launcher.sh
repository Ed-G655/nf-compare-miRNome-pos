##!/usr/bin/env bash
cd ../
ref_dir="/bodega/projects/miRNome_project/analisis/100GMX/targets/ref_targets/"
alt_dir="/bodega/projects/miRNome_project/analisis/100GMX/targets/alt_targets/"
bed="/bodega/projects/miRNome_project/analisis/100GMX/hsa.bed"
output_directory="$(dirname $ref_dir)/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run nf-compare-miRNome-pos.nf \
	--ref_dir $ref_dir \
  --alt_dir $alt_dir \
	--bed $bed \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
	-with-timeline $output_directory/`date +%Y%m%d_%H%M%S`_timeline.html \
  && echo -e "======\n  pipeline execution END \n======"
