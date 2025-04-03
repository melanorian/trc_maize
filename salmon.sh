#!/bin/bash

# Replace with the actual paths to your directories and index
output_dir=/home/melanie/net/virus_melanie/maize/results/t4_all_v1/  # Output directory for Salmon quantification
R1=/home/melanie/net/virus_melanie/maize/raw_data/demultiplexed/R1  # Path to R1 reads
R2=/home/melanie/net/virus_melanie/maize/raw_data/demultiplexed/R2  # Path to R2 reads
index=/home/melanie/net/virus_melanie/maize/annotation/salmon_index/index  # Salmon index directory
threads=50  # Number of threads
log_file=${output_dir}/salmon_quant.log  # Log file for tracking

# Create output directory if it doesn't exist
if [ ! -d ${output_dir} ]; then
    mkdir -p ${output_dir}
    log_message "Output directory created: ${output_dir}"
fi

cd "${output_dir}"

# Function to log messages
log_message() {
    echo $(date +'%Y-%m-%d %H:%M:%S') $1 >> ${log_file}
}

# Check if directories exist
if [ ! -d ${R1} ]; then
    log_message "Error: R1 does not exist."
    exit 1
fi

if [ ! -d ${R2} ]; then
    log_message "Error: R2 does not exist."
    exit 1
fi

if [ ! -d ${index} ]; then
    log_message "Error: Salmon index does not exist."
    exit 1
fi

log_message "Directories and files exist. Starting processing..."

# Make a file listing all the individual sample names (before the last underscore)
ls ${R1} | awk -F'_' '{print $1}' | sort -u > ${output_dir}/samplenames.txt

# Read each line from samplenames.txt and process the corresponding files
cat ${output_dir}/samplenames.txt | while read -r line
do
    # Extract the file paths for R1 and R2 based on the line (sample name)
    file1=${R1}/${line}_1.fq.gz  # R1 paired-end file
    file2=${R2}/${line}_2.fq.gz  # R2 paired-end file

    # Check if files exist for the current sample
    if [ ! -f ${file1} ] || [ ! -f ${file2} ]; then
        log_message "Error: Files not found for sample ${line}. Skipping..."
        continue
    fi

    log_message "Processing ${file1} ${file2}"  

    # Run Salmon quantification (mapping and quantification in one step)
    /home/melanie/tools/anaconda3/envs/salmon/bin/salmon quant \
        -i ${index} \
        -l A \
        -1 ${file1} \
        -2 ${file2} \
        --validateMappings \
        -o ${output_dir}/${line} \
        -p ${threads} \
        >> ${log_file} 2>&1
done

log_message "Salmon quantification completed for all samples."