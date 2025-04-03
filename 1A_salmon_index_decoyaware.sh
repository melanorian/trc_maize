#!/bin/bash

# user input
output_dir=/home/melanie/net/virus_melanie/maize/annotation/salmon_index/
gentrome=/home/melanie/net/virus_melanie/maize/assembly/Zm-B73-gentrome.fa
decoys=/home/melanie/net/virus_melanie/maize/assembly/decoys.txt
threads=50

cd ${output_dir}

# Check if input files exist
if [ ! -f ${gentrome} ] || [ ! -f ${decoys} ]; then
  echo "Error: One or more input files not found."
  exit 1
fi

# run Salmon index function
salmon index -t ${gentrome} -d ${decoys} -k 31 -p ${threads} -i ${output_dir}/index

cd
