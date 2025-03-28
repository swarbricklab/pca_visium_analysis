#!/bin/bash

set -euxo pipefail

# --- Aim: submit jobs for a certain script for all samples (on wolfpack hpc)

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

run_script_file="$cwd/run_universal.sh" # downstream bash script that will run the python script
function_script_file="$cwd/$file_name" # actual python script to be run


# --- Identify samples
# Identify all from the sample sheet
# Extract the section_id column
samples=$(awk -F',' 'NR > 1 {print $5}' $visium_sample_sheet) 

# samples=("3946_2P_Xe1_XV" "3962_3H_Xe1_XV" "4497-1_2K_Xe1_XV" "4754_1G_Xe1_XV")
# samples=("3946_2P_Xe1_V" "3962_3H_Xe1_V" "4497-1_2K_Xe1_V" "4754_1G_Xe1_V")

echo "Identified samples: $samples"

# --- Submit jobs

n_cores=4

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
        -M s.vanderleij@garvan.org.au \
        -m a \
        -V \
        -pe smp $n_cores \
        -l mem_requested=8G \
        $run_script_file $function_script_file $raw_data_dir $h5ad_dir $project_dir $s $environment $min_gene $min_count
        
done
