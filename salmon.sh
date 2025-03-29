#!/bin/bash

# STAR + Salmon pipeline for multiple RNA-Seq samples with decoy-aware quantification

# User input paths
genome_fasta=/home/melanie/net/virus_melanie/maize/assembly/Zm-B73-REFERENCE-NAM-5.0.fa  # Genome fasta file (uncompressed)
genome_gff=/home/melanie/net/virus_melanie/maize/annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3  # GFF annotation file
genome_dir=/home/melanie/net/virus_melanie/maize/annotation/STAR_index  # STAR genome index directory
output_dir=/home/melanie/net/virus_melanie/maize/results/t2_star_salmon  # STAR output directory
salmon_index=/home/melanie/net/virus_melanie/maize/annotation/salmon_index  # Salmon index directory
gentrome_fasta=/home/melanie/net/virus_melanie/maize/annotation/gentrome.fa  # Gentrome (transcriptome + genome concatenated)
decoys_file=/home/melanie/net/virus_melanie/maize/annotation/decoys.txt  # Decoy sequences file (chromosome names)
transcriptome_fasta=/home/melanie/net/virus_melanie/maize/annotation/transcripts.fa  # Transcriptome for Salmon

# Software paths
STAR=/home/melanie/tools/STAR/STAR-2.7.11a/source/STAR  # STAR executable
log_file=${output_dir}/log_file.txt  # Log file

# Create output directories if they don't exist
mkdir -p ${output_dir}
mkdir -p ${salmon_index}

# Create log file if it doesn't exist
if [ ! -f ${log_file} ]; then
    touch ${log_file}
fi

# Set threads
threads=50

# Check if STAR genome index exists
if [ ! -d ${genome_dir} ]; then
    echo "$(date): Genome index not found. Generating genome index..." >> ${log_file}
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

# Check if Salmon index exists
if [ ! -d ${salmon_index}/ ]; then
    echo "$(date): Salmon index not found. Generating Salmon index..." >> ${log_file}
    
    # Generate Salmon index using the gentrome (genome + transcriptome concatenated)
    salmon index -t ${gentrome_fasta} -d ${decoys_file} -i ${salmon_index}

    echo "$(date): Salmon index creation completed." >> ${log_file}
else
    echo "$(date): Salmon index already exists. Skipping index generation." >> ${log_file}
fi

# Loop over each sample in the R1 directory
for sample in /home/melanie/net/virus_melanie/maize/raw_data/demultiplexed/R1/*_1.fq.gz; do
    # Extract sample name
    sample_name=$(basename ${sample} "_1.fq.gz")
    sample_name=$(echo ${sample_name} | sed 's/^_//;s/_$//')

    # Define file paths
    rnaseq_R1=/home/melanie/net/virus_melanie/maize/raw_data/demultiplexed/R1/${sample_name}_1.fq.gz
    rnaseq_R2=/home/melanie/net/virus_melanie/maize/raw_data/demultiplexed/R2/${sample_name}_2.fq.gz
    bam_file=${output_dir}/${sample_name}_Aligned.sortedByCoord.out.bam
    salmon_output=${output_dir}/salmon_quant/${sample_name}

    # Check if R2 file exists
    if [ -f ${rnaseq_R2} ]; then
        echo "$(date): Processing sample ${sample_name}" >> ${log_file}

        # Check if BAM file already exists before running STAR
        if [ ! -f ${bam_file} ]; then
            # Run STAR for alignment if BAM file does not exist
            ${STAR} --runMode alignReads --genomeDir ${genome_dir} \
            --readFilesIn ${rnaseq_R1} ${rnaseq_R2} --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
            --runThreadN ${threads} \
            --outFileNamePrefix ${output_dir}/${sample_name}_

            echo "$(date): STAR mapping completed for ${sample_name}" >> ${log_file}
        else
            echo "$(date): BAM file already exists for ${sample_name}. Skipping STAR alignment." >> ${log_file}
        fi

        # Check if quant.sf file already exists for the sample
        if [ ! -f ${salmon_output}/quant.sf ]; then
            # Run Salmon for quantification (decoy-aware mode)
            mkdir -p ${salmon_output}
            /home/melanie/tools/anaconda3/envs/salmon/bin/salmon quant -i ${salmon_index} -l A -p ${threads} --gcBias -a ${bam_file} -o ${salmon_output}

            echo "$(date): Salmon quantification completed for ${sample_name}" >> ${log_file}
        else
            echo "$(date): quant.sf already exists for ${sample_name}. Skipping Salmon quantification." >> ${log_file}
        fi

    else
        echo "$(date): R2 file for ${sample_name} not found. Skipping." >> ${log_file}
    fi
done

echo "$(date): STAR + Salmon pipeline completed for all samples." >> ${log_file}
