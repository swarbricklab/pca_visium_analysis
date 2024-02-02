#!/bin/bash

# --- Aim: automatically a python or R script on multiple samples

# print the date & time for the log file
now=$(date)
echo
echo "Start time: $now"

# arguments from submission file
function_script_file=$1
h5ad_dir=$2
project_dir=$3
environment=$4
res_dir=$5
use_data=$6
shift 6 # need to shift so that the first 6 arguments are gone and I get only the array of samples now
samples=$@

# load environment
conda activate $environment
echo "Run script $function_script_file"

# script file 
$function_script_file \
    --h5ad_dir $h5ad_dir \
    --project_dir $project_dir \
    --res_dir $res_dir \
    --use_data $use_data \
    --samples $samples \
