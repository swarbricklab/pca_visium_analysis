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
# source "../config/05_train_model_major_params.sh"
# source "../config/08_train_model_minor_mal_params.sh"
# source "../config/16_train_model_major_temp_params.sh"

# source "../config/03a_celltype_mapping_params.sh"
# source "../config/06a_celltype_mapping_major_params.sh"
# source "../config/09a_celltype_mapping_minor_mal_params.sh"
# source "../config/17a_celltype_mapping_major_temp_params.sh"

# source "../config/04_visualize_results_params.sh"
source "../config/07_visualize_results_major_params.sh"
# source "../config/14_visualize_results_minor_mal_params.sh"
# source "../config/18_visualize_results_major_temp_params.sh"

logDir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/02_cell2location/logs/"
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
samples=$(awk -F',' 'NR > 1 {print $2}' $visium_sample_sheet)

# samples=('20111-1' '20384-2' '20216-1' '19617-2' '20153-1' '20153-2' '20033' '20111-2' '20130-1' '20130-2', 'PCa20153_C1_20128_C1', 'PCa20130_C1_20272_C2')
# samples=('PCa20130_C1_20272_C2')
# samples='C2L_major_temp' # if running as a single job
# samples=( "H1_2" "H1_4" "H1_5" "H2_1" "H2_2" "H2_5" "V1_2")

# submit a single job for training the model

# --- Submit jobs

n_cores=4
# n_cores=8

# Use script name as name for job
s=$file_name
tmp=${s##*/}   # trim the directory path (split on the last '/' and take the last element)
tmp=${tmp##*_}   # also trim the digit at the start
name=${tmp%%.*}  # get rid of .py (split on '.' and take first element)

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
        $run_script_file $function_script_file $projectName $repo $exp $analysis $s $environment $h5ad_dir $model_dir
done

### FOR TRAINING MODEL: $run_script_file $function_script_file $projectName $repo $exp $analysis $s $environment
### FOR RUNNING C2L: $run_script_file $function_script_file $projectName $repo $exp $analysis $s $environment $h5ad_dir $model_dir