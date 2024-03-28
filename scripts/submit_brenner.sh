#!/bin/bash

source /share/ScratchGeneral/sopvan/conda/envs/activate cell2location

set -euxo pipefail

# --- Aim: submit jobs that need to use a gpu (on Brenner)

while getopts s:p: flag
do
    case "${flag}" in
        s) file_name=${OPTARG};; # name of the script you want to run (inlcuding the path to it from the current directory)
    esac
done

echo "Run script $file_name"

# --- Load parameters

# source "../config/02_train_model_params.sh"
# source "../config/03a_celltype_mapping_params.sh"
# source "../config/04_visualize_results_params.sh"
source "../config/05_train_model_major_params.sh"

logDir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/03_mengxiao_data_PMID_35948708/logs/"
mkdir -p $logDir

# --- Paths to donwstream scripts
cwd=$(pwd)

run_script_file="$cwd/run_universal.sh" # downstream bash script that will run the python script
function_script_file="$cwd/$file_name" # actual python script to be run


# --- Identify samples
# Identify all from the sample sheet
# extract the unique sample id from the sample sheet
# Extract the section_id column
# samples=$(awk -F',' 'NR > 1 {print $5}' $visium_sample_sheet) 

# echo "Identified samples: $samples"

# submit a single job for training the model

# --- Submit jobs

n_cores=8

# Use script name as name for job
s=$file_name
tmp=${s##*/}   # trim the directory path (split on the last '/' and take the last element)
tmp=${tmp##*_}   # also trim the digit at the start
name=${tmp%%.*}  # get rid of .py (split on '.' and take first element)

# samples=('20111-1' '20384-2' '20216-1' '19617-2' '20153-1' '20153-2' '20033' '20111-2' '20130-1' '20130-2')
samples='C2L_major_all' # if running as a single job


for s in ${samples[@]}
do 
    #cd $repoDir
    echo $s

    qsub -N "${name}_${s}" \
        -cwd \
        -S /bin/bash \
        -o $logDir/$s.out \
        -e $logDir/$s.err \
        -M e.apostolov@garvan.org.au \
        -m a \
        -V \
        -pe smp $n_cores \
        -l mem_requested=12G \
        -l nvgpu=1 \
        -l h='epsilon*' \
        $run_script_file $function_script_file $projectName $repo $exp $analysis $s $environment

done

###        $run_script_file $function_script_file $projectName $repo $exp $analysis $s $environment $h5ad_dir $model_dir
