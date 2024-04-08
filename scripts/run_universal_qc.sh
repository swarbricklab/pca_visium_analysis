#!/bin/bash

# --- Aim: automatically a python or R script on a specific sample

# print the date & time for the log file
now=$(date)
echo
echo "Start time: $now"

# arguments from submission file
function_script_file=$1
h5ad_dir=$2
project_dir=$3
sample_id=$4
environment=$5
min_count=$6
min_gene=$7
min_spots=$8
res_dir=$9

# load environment
conda activate $environment
echo "Run script $function_script_file"

# script file 
$function_script_file \
    --h5ad_dir $h5ad_dir \
    --sample_id $sample_id \
    --project_dir $project_dir \
    --min_count $min_count \
    --min_gene $min_gene \
    --min_spots $min_spots \
    --res_dir $res_dir \