#!/bin/bash

# --- Aim: automatically a python or R script on a specific sample

# print the date & time for the log file
now=$(date)
echo
echo "Start time: $now"

# arguments from submission file
function_script_file=$1
raw_data_dir=$2
h5ad_dir=$3
project_dir=$4
sample_id=$5
environment=$6
min_gene=$7
min_count=$8

# load environment
conda activate $environment
echo "Run script $function_script_file"

# script file 
$function_script_file \
    --raw_data_dir $raw_data_dir \
    --h5ad_dir $h5ad_dir \
    --sample_id $sample_id \
    --project_dir $project_dir \
    --min_gene $min_gene \
    --min_count $min_count \
