#!/bin/bash

#SBATCH --job-name=Repeatome
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jlamb@crimson.ua.edu
#SBATCH -N 1 #nodes
#SBATCH -o Repeatome_%j.o

echo "Starting the pipeline at $(date)"
#usage = Repeatome_Pipeline3.sh <fasta_file>
# The input fasta file should be paired end and interleaved with duplicates removed

Fasta=$1
Base=$(basename "$Fasta" .fasta)

echo "Processing samples with base name: $Base"

ProjectDir="/home/jlamb1/Projects/Repeatome/$Base"

[ ! -d "$ProjectDir" ] && mkdir -p "$ProjectDir"

cp "$Fasta" "$ProjectDir"

cd "$ProjectDir" || { echo "Failed to change directory to $ProjectDir"; exit 1; }
#copy the fasta file to the project directory

zero_SATS_count=0
zero_LTRS_count=0
Num_of_runs=0

while :; do
    Num_of_runs=$((Num_of_runs + 1))

    # Take a random sample of 1 million reads
    echo "Sampling 1 million reads..."
    seqtk sample -s100 ${Base}.fasta 1000000 > ${Base}_sample.fasta

    # Run repeat explorer
    echo "Running Repeat Explorer..."
    source activate /home/jlamb1/bin/miniconda3/envs/eccsplorer
    /home/jlamb1/bin/repex_tarean/seqclust -p -A -v ./re_output -s 1000000 -c 32 -C -tax METAZOA3.0 -D BLASTX_W2 ${Base}_sample.fasta

    # Check if Repeat Explorer was successful
    if [ $? -ne 0 ]; then
        echo "Error during Repeat Explorer execution. Exiting."
        exit 1
    fi

    # Deactivate the conda environment
    conda deactivate

    # Find the number of Repeats found
    high_conf_SATs=$(grep -c "^>" ./re_output/TAREAN_consensus_rank_1.fasta)
    high_conf_LTRs=$(grep -c "^>" ./re_output/TAREAN_consensus_rank_3.fasta)

    echo "There were $high_conf_SATs high confidence SATs found"
    echo "There were $high_conf_LTRs high confidence LTRs found"

    # Check if num of repeats found is > 0
    if [ $high_conf_SATs -gt 0 ]; then
        zero_SATS_count=0  # reset count if non-zero found
    else
        zero_SATS_count=$((zero_SATS_count + 1))
    fi
    if [ $high_conf_LTRs -gt 0 ]; then
        zero_LTRS_count=0  # reset count if non-zero found
    else
        zero_LTRS_count=$((zero_LTRS_count + 1))
    fi

    # Check if two consecutive iterations returned 0 for both SATS and LTRS or if the loop has run for more than, say, 10 times
    if [ $zero_SATS_count -ge 2 ] && [ $zero_LTRS_count -ge 2 ] || [ $Num_of_runs -gt 10 ]; then
        echo "Two consecutive iterations with zero high confidence SATs and LTRs found, or max run limit reached. Exiting."
        break
    fi

    # Check if a repeat_library_fasta file exists, if not make one
    LibraryFile="$ProjectDir/repeat_library.fasta"
    if [ ! -f "$LibraryFile" ]; then
        touch "$LibraryFile"
    fi

    # Activate conda env
    source activate /home/jlamb1/bin/miniconda3/envs/bowtie2

    # Append the files to the library
    cat ./re_output/TAREAN_consensus_rank_1.fasta ./re_output/TAREAN_consensus_rank_3.fasta >> "$LibraryFile"

    # Index the reads
    bwa index ${Base}.fasta
    echo "Mapping reads to the library..."

    # Map the reads
    bwa mem ${Base}.fasta "$LibraryFile" > ${Base}_${Num_of_runs}_alignment.sam
    echo "Done mapping reads to the library"

    # Convert SAM to BAM
    samtools view -S -b ${Base}_${Num_of_runs}_alignment.sam > ${Base}_${Num_of_runs}_alignment.bam
    echo "Done converting SAM to BAM"

    # Sort the BAM file
    samtools sort ${Base}_${Num_of_runs}_alignment.bam -o ${Base}_${Num_of_runs}_sorted_alignment.bam
    echo "Done sorting BAM file"

    # Filter out mapped reads
    samtools view -b -F 4 ${Base}_${Num_of_runs}_sorted_alignment.bam > ${Base}_${Num_of_runs}_mapped_reads.bam
    echo "Done filtering out mapped reads"

    # Get list of mapped read names
    samtools view ${Base}_${Num_of_runs}_mapped_reads.bam | cut -f1 > ${Base}_${Num_of_runs}_mapped_read_names.txt
    echo "Done getting list of mapped read names"
    # Use seqtk to filter out mapped reads from the original FASTA
    seqtk subseq ${Base}.fasta ${Base}_${Num_of_runs}_mapped_read_names.txt > ${Base}_${Num_of_runs}_filtered.fasta
    echo "Done filtering out mapped reads from the original FASTA"

    # Deactivate conda
    conda deactivate


    mv ${Base}_${Num_of_runs}_filtered.fasta ${Base}.fasta 
    echo "Done moving filtered FASTA to original FASTA"
    
    rm ${Base}_${Num_of_runs}_alignment.sam ${Base}_${Num_of_runs}_alignment.bam ${Base}_${Num_of_runs}_sorted_alignment.bam ${Base}_${Num_of_runs}_mapped_reads.bam ${Base}_${Num_of_runs}_mapped_read_names.txt

    echo "Temp files removed"
done
