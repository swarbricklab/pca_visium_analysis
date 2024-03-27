#!/bin/bash
R="/home/evaapo/miniconda3/envs/convert_data/bin/R"

# samples
samples=( "H1_2" "H1_4" "H1_5" "H2_1" "H2_2" "H2_5" "V1_2")

for sample in ${samples[@]}; do

	email="e.apostolov@garvan.org.au"
	cores=8
	group="TumourProgression"

	# set up directory structure
	projectName="PCa_Visium"
	# repo name
	repo="pca_visium_analysis"
	# experiment code
	exp="03_mengxiao_data_PMID_35948708"
	# analysis
	analysis="01_create_h5ad_objects"
	# samples
	sample_id=$sample

	#seuratScript=$projectDir"/scripts/atac/${analysis}/${run}.R"
	seuratScript="/share/ScratchGeneral/evaapo/projects/${projectName}/${repo}/scripts/${exp}/${analysis}.R"

	jobName="${sample}_${analysis}"
	logDir="/share/ScratchGeneral/evaapo/projects/${projectName}/${repo}/results/${exp}/${analysis}/logs/${jobName}/"
	#paramFile="/share/ScratchGeneral/evaapo/projects/${projectName}/scripts/atac/${analysis}/01_params_file_2.csv"
	mkdir -p $logDir
	qsub -wd $logDir -q short.q -pe smp $cores -l mem_requested=8G -l tmp_requested=200G -m ea -M $email -b y -j y -V -P $group -N $jobName "${R} CMD BATCH --no-save '--args $projectName $exp $repo $analysis $sample_id' $seuratScript"

done