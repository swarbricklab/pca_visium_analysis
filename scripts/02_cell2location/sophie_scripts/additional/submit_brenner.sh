#!/bin/bash

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

source "../config/params.sh"


# --- Paths to donwstream scripts
cwd=$(pwd)

run_script_file="$cwd/run_universal.sh" # downstream bash script that will run the python script
function_script_file="$cwd/$file_name" # actual python script to be run


# --- Identify samples
# Identify all from the sample sheet
# extract the unique sample id from the sample sheet
# Extract the section_id column
samples=$(awk -F',' 'NR > 1 {print $5}' $visium_sample_sheet) 

echo "Identified samples: $samples"

# submit a single job for training the model
# samples='cell2location_model'
# samples='4754_1G_Xe1_XV'
# samples=('3946_2P_Xe1_V' '3946_2P_Xe1_XV' '3962_3H_Xe1_V')
samples='cell2loc_all'


# --- Submit jobs

n_cores=8

# Use script name as name for job

s=$file_name
tmp=${s##*/}   # trim the directory path (split on the last '/' and take the last element)
tmp=${tmp##*_}   # also trim the digit at the start
name=${tmp%%.*}  # get rid of .py (split on '.' and take first element)

# Set environment to cell2location
environment='cell2location'


for s in ${samples[@]}
do 
    cd $project_dir
    echo $s 

    qsub -N "${name}_${s}" \
        -cwd \
        -S /bin/bash \
        -o $log_dir/$s.out \
        -e $log_dir/$s.err \
        -M s.vanderleij@garvan.org.au \
        -m a \
        -V \
        -pe smp $n_cores \
        -l mem_requested=12G \
        -l nvgpu=1 \
        -l h='epsilon*' \
        $run_script_file $function_script_file $raw_data_dir $h5ad_dir $project_dir $s $environment $min_gene $min_count
        
done


