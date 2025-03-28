#!/bin/bash

set -euxo pipefail

# --- Aim: submit jobs for a certain script for all samples

while getopts s:p: flag
do
    case "${flag}" in
        s) file_name=${OPTARG};; # name of the script you want to run (inlcuding the path to it from the current directory)
    esac
done

echo "Run script $file_name"

# --- Load parameters

# source "../config/params_qc.sh"
source "../config/params_qc.sh"

# --- Paths to donwstream scripts

run_script_file="$project_dir/scripts/run_universal_qc.sh" # downstream bash script that will run the python script
function_script_file="$project_dir/scripts/$file_name" # actual python script to be run

# --- Identify samples
# Identify all from the sample sheet
# Extract the section_name column
samples=$(awk -F',' 'NR > 1 {print $2}' $visium_sample_sheet)
# samples='All_samples' # if running as a single job

echo "Identified samples: $samples"

# --- Submit jobs

n_cores=8

# Use script name as name for job

s=$file_name
tmp=${s##*/}   # trim the directory path (split on the last '/' and take the last element)
tmp=${tmp##*_}   # also trim the digit at the start
name=${tmp%%.*}  # get rid of .py (split on '.' and take first element)

for s in ${samples[@]}
do 
    cd $project_dir
    echo $s 

    qsub -N "${name}_${s}" \
        -cwd \
        -S /bin/bash \
        -o $log_dir/$s.out \
        -e $log_dir/$s.err \
        -M e.apostolov@garvan.org.au \
        -m a \
        -V \
        -pe smp $n_cores \
        -l mem_requested=8G \
        $run_script_file $function_script_file $h5ad_dir $project_dir $s $environment $min_count $min_gene $min_spots $res_dir

done

# for all plots QC - doesn't work, issue with environemnt: $run_script_file $function_script_file $h5ad_dir $project_dir $s $environment $res_dir
