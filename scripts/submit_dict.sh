#!/bin/bash

set -euxo pipefail

# --- Aim: submit a job for all visium samples + a particular script

while getopts s:p: flag
do
    case "${flag}" in
        s) file_name=${OPTARG};; # name of the script you want to run (inlcuding the path to it from the current directory)
    esac
done

echo "Run script $file_name"

# --- Load parameters

source "../config/params.sh"


# --- Paths to donwstream scripts
run_script_file="$cwd/run_dict.sh" # downstream bash script that will run the python script
function_script_file="$cwd/$file_name" # actual python script to be run


# --- Identify samples
# Identify all from the sample sheet
# extract the unique sample id from the sample sheet
samples=$(awk -F',' 'NR > 1 {print $3}' $visium_sample_sheet) 

echo "Identified samples: $samples"

# --- Submit jobs
log_dir="$project_dir/logs"
mkdir -p $log_dir

n_cores=8

# Use script name as name for job

s=$file_name
tmp=${s##*/}   # trim the directory path (split on the last '/' and take the last element)
tmp=${tmp##*_}   # also trim the digit at the start
name=${tmp%%.*}  # get rid of .py (split on '.' and take first element)

log_name=$name

cd $project_dir
echo $samples

qsub -N "${name}_all_samples" \
    -cwd \
    -S /bin/bash \
    -o $log_dir/$log_name.out \
    -e $log_dir/$log_name.err \
    -M s.vanderleij@garvan.org.au \
    -m a \
    -V \
    -pe smp $n_cores \
    -l mem_requested=8G \
    $run_script_file $function_script_file $h5ad_dir $project_dir $environment $res_dir $use_data ${samples[@]}


