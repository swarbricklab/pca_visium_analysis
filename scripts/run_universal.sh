#!/bin/bash


# --- Aim: automatically a python or R script on a specific sample

# print the date & time for the log file
now=$(date)
echo
echo "Start time: $now"

# arguments from submission file
function_script_file=$1
projectName=$2
repo=$3
exp=$4
analysis=$5
sample_id=$6
environment=$7

# Paths to downstream scripts
h5ad_dir=$8
model_dir=$9

# load environment
conda activate "$environment"
echo "Run script $function_script_file"

# script file 
$function_script_file \
    --projectName $projectName \
    --repo $repo \
    --exp $exp \
    --analysis $analysis \
    --sample_id $sample_id \
    --h5ad_dir $h5ad_dir \
    --model_dir $model_dir \