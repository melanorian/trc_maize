#!/bin/bash

# STAR run for multiple demultiplexed samples (R1 and R2)

# user input
genome_fasta=/home/melanie/net/virus_melanie/maize/assembly/Zm-B73-REFERENCE-NAM-5.0.fa  # Genome fasta file (uncompressed)
genome_gff=/home/melanie/net/virus_melanie/maize/annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3  # GFF annotation file
genome_dir=/home/melanie/net/virus_melanie/maize/annotation/STAR_index  # STAR genome index directory
output_dir=/home/melanie/net/virus_melanie/maize/results/star_alignment  # Output directory
STAR=/home/melanie/tools/STAR/STAR-2.7.11a/source/STAR  # STAR executable
log_file=/home/melanie/net/virus_melanie/maize/results/star_alignment/log_file.txt  # Log file to track processing

# Create output directory if it doesn't exist
if [ ! -d ${output_dir} ]; then
    mkdir -p ${output_dir}
fi

# Create log file if it doesn't exist
if [ ! -f ${log_file} ]; then
    touch ${log_file}
fi

# Set threads and maximum intron size
threads=50
max_intron=900

# Check if genome index already exists
if [ ! -d ${genome_dir} ]; then
    echo "$(date): Genome index not found. Generating genome index..." >> ${log_file}

    # Create the STAR genome index directory
    mkdir -p ${genome_dir}

    # Index genome assembly with STAR
    ${STAR} --runMode genomeGenerate --genomeDir ${genome_dir} \
    --genomeFastaFiles ${genome_fasta} --sjdbGTFfile ${genome_gff} \
    --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene Parent --sjdbGTFfeatureExon CDS \
    --runThreadN ${threads} --genomeSAindexNbases 11

    echo "$(date): Genome index generation completed." >> ${log_file}
else
    echo "$(date): Genome index already exists. Skipping index generation." >> ${log_file}
fi

# Loop over each sample in the R1 directory
for sample in /home/melanie/net/virus_melanie/maize/raw_data/demultiplexed/R1/*_1.fq.gz; do
    # Extract sample name from the file (removing the _1.fq.gz part)
    sample_name=$(basename ${sample} "_1.fq.gz")

    # Remove trailing or leading underscores from sample name (if any)
    sample_name=$(echo ${sample_name} | sed 's/^_//;s/_$//')

    # Set the corresponding R1 and R2 file paths
    rnaseq_R1=/home/melanie/net/virus_melanie/maize/raw_data/demultiplexed/R1/${sample_name}_1.fq.gz
    rnaseq_R2=/home/melanie/net/virus_melanie/maize/raw_data/demultiplexed/R2/${sample_name}_2.fq.gz

    # Check if the corresponding R2 file exists
    if [ -f ${rnaseq_R2} ]; then
        echo "$(date): Processing sample ${sample_name}" >> ${log_file}

        # Run STAR for alignment and quantification
        ${STAR} --runMode alignReads --genomeDir ${genome_dir} \
        --readFilesIn ${rnaseq_R1} ${rnaseq_R2} --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
        --quantMode GeneCounts --runThreadN ${threads} \
        --outFileNamePrefix ${output_dir}/${sample_name}_  # this line

        echo "$(date): Finished processing sample ${sample_name}" >> ${log_file}

    else
        echo "$(date): R2 file for sample ${sample_name} not found. Skipping this sample." >> ${log_file}
    fi
done

# Completion message
echo "$(date): STAR alignment completed for all samples." >> ${log_file}
