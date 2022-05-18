ref_dir="test/data/ref_targets/"
alt_dir="test/data/alt_targets/"
bed="test/data/sample.bed"
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
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
	#-stub-run \
