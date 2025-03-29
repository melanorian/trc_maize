#!/bin/bash

# user input
output_dir=/net/virus/linuxhome/michael-group/melanie/transcriptomics_timeseries/results/salmon_combined/
gentrome=/net/virus/linuxhome/michael-group/melanie/transcriptomics_timeseries/raw_data/Salmon_input/MM20231127_Pe1_Sp_gentrom_v2.fa
decoys=/net/virus/linuxhome/michael-group/melanie/transcriptomics_timeseries/raw_data/Salmon_input/MM20231127_Pe_Sp_decoy.txt
threads=19

cd ${output_dir}

# Check if input files exist
if [ ! -f ${gentrome} ] || [ ! -f ${decoys} ]; then
  echo "Error: One or more input files not found."
  exit 1
fi

# run Salmon index function
salmon index -t ${gentrome} -d ${decoys} -k 31 -p ${threads} -i ${output_dir}/index

cd
