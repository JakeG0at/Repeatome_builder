#!/bin/bash
#SBATCH --job-name=Repeatome
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jlamb@crimson.ua.edu
#SBATCH -N 1 #nodes
#SBATCH -o Repeatome_%j.o
threads=32
#divide the number of threads by 2
threads2=$((threads / 2))

echo "Starting the pipeline at $(date)"
# I will hard code the file names for now

#copy the raw fastq files to the project directory
#cp /home/jlamb1/SRA_seqs/SRX19958881_N1.fastq.gz /home/jlamb1/Projects/Repeatome/SRX19958881/ &
#cp /home/jlamb1/SRA_seqs/SRX19958881_N2.fastq.gz /home/jlamb1/Projects/Repeatome/SRX19958881/ &
#wait

#gunzip the files if they are gzipped
#if [ -f /home/jlamb1/Projects/Repeatome/SRX19958881/SRX19958881_N1.fastq.gz ]; then
   #echo "Unzipping the files"
    #gunzip /home/jlamb1/Projects/Repeatome/SRX19958881/SRX19958881_N1.fastq.gz &
   # gunzip /home/jlamb1/Projects/Repeatome/SRX19958881/SRX19958881_N2.fastq.gz &
    #wait
  #  echo "Done unzipping the files"
#fi

read1="/home/jlamb1/Projects/Repeatome/SRX19958881/SRX19958881_N1.fastq"
read2="/home/jlamb1/Projects/Repeatome/SRX19958881/SRX19958881_N2.fastq"

#take the basename of the file and use that as the base name for the project directory
Base=$(basename "$read1" _N1.fastq)
echo "Processing samples with base name: $Base"

ProjectDir="/home/jlamb1/Projects/Repeatome/$Base"
[ ! -d "$ProjectDir" ] && mkdir -p "$ProjectDir"
cd "$ProjectDir" || { echo "Failed to change directory to $ProjectDir"; exit 1; }
mkdir -p "$ProjectDir"/Repeat_explorer_outputs
zero_SATS_count=0
zero_LTRS_count=0
Num_of_runs=0

while :; do
    Num_of_runs=$((Num_of_runs + 1))

    echo "Fixing the reads..."
    /home/jlamb1/bin/bbmap/repair.sh in1=$read1 in2=$read2 out1=${Base}_N1_fixed.fastq out2=${Base}_N2_fixed.fastq outs=${Base}_singletons.fastq ow=true threads=$threads
    echo "Done fixing the reads"

    # Update the read1 and read2 variables to point to the fixed files
    read1="${Base}_N1_fixed.fastq"
    read2="${Base}_N2_fixed.fastq"

    # Take a random sample of 1 million reads
    echo "Sampling 1 million reads..."
    seqtk sample -s100 $read1 500000 > sub1.fastq &
    seqtk sample -s100 $read2 500000 > sub2.fastq &
    wait

    # Get the headers from the sampled reads
    grep "^@" sub1.fastq | sed 's/^@//' > headers1.txt &
    grep "^@" sub2.fastq | sed 's/^@//' > headers2.txt &
    wait
    # Combine the headers from both files
    cat headers1.txt headers2.txt > headers.txt

    #remove matching headers from read1 and read2
    grep -v -F -f headers.txt $read1 > filtered1.fastq &
    grep -v -F -f headers.txt $read2 > filtered2.fastq &
    wait
    #Change file names to original names
    mv filtered1.fastq $read1
    mv filtered2.fastq $read2

    #Convert to fasta that is merged and interleaved
    echo "Converting to interleaved fasta"
    seqtk mergepe sub1.fastq sub2.fastq > merged.fastq
    seqtk seq -A merged.fastq > ${Base}_sample.fasta

    # Run repeat explorer
    echo "Running Repeat Explorer..."
    source activate /home/jlamb1/bin/miniconda3/envs/eccsplorer || { echo "Failed to activate conda environment"; exit 1; }

    /home/jlamb1/bin/repex_tarean/seqclust -p -A -v ./re_output -s 1000000 -c $threads -C -tax METAZOA3.0 -D BLASTX_W2 ${Base}_sample.fasta
    # Check if Repeat Explorer was successful
    if [ $? -eq 0 ]; then
        echo "Repeat Explorer was successful"
    else
        echo "Repeat Explorer failed"
        exit 1
    fi

    # Deactivate the conda environment
    conda deactivate || { echo "Failed to deactivate conda environment"; exit 1; }

    # Find the number of Repeats found
    high_conf_SATs=$(grep -c "^>" ./re_output/TAREAN_consensus_rank_1.fasta) &
    high_conf_LTRs=$(grep -c "^>" ./re_output/TAREAN_consensus_rank_3.fasta) &
    wait

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

    # Check if two consecutive iterations returned 0 for both SATS and LTRS
    if [ $zero_SATS_count -ge 2 ] && [ $zero_LTRS_count -ge 2 ]; then
        echo "Two consecutive iterations with zero high confidence SATs and LTRs found, or max run limit reached. Exiting."
        break
    fi

    # Check if a repeat_library_fasta file exists, if not make one
    LibraryFile="$ProjectDir/repeat_library.fasta"
    if [ ! -f "$LibraryFile" ]; then
        touch "$LibraryFile"
    fi

    # Activate conda env
    source activate /home/jlamb1/bin/miniconda3/envs/bowtie2 || { echo "Failed to activate conda environment"; exit 1; }
    # Append the files to the library
    cat ./re_output/TAREAN_consensus_rank_1.fasta ./re_output/TAREAN_consensus_rank_3.fasta >> "$LibraryFile"

    # Ensure the reference library is indexed
echo "Indexing the reference library..."
bwa index "$LibraryFile"
echo "Done indexing the reference library"

echo "Mapping reads to the library..."
# Map read1 and read2 to the library
bwa mem -t $threads2 "$LibraryFile" $read1 > ${Base}_${Num_of_runs}_alignment1.sam &
bwa mem -t $threads2 "$LibraryFile" $read2 > ${Base}_${Num_of_runs}_alignment2.sam &
wait
echo "Done mapping reads to the library"

# Convert SAM to BAM
echo "Converting SAM to BAM..."
samtools view -@ $threads2 -S -b ${Base}_${Num_of_runs}_alignment1.sam > ${Base}_${Num_of_runs}_alignment1.bam &
samtools view -@ $threads2 -S -b ${Base}_${Num_of_runs}_alignment2.sam > ${Base}_${Num_of_runs}_alignment2.bam &
wait
echo "Done converting SAM to BAM"

# Sort the BAM file
echo "Sorting BAM file..."
samtools sort -@ $threads2 ${Base}_${Num_of_runs}_alignment1.bam -o ${Base}_${Num_of_runs}_alignment1_sorted.bam &
samtools sort -@ $threads2 ${Base}_${Num_of_runs}_alignment2.bam -o ${Base}_${Num_of_runs}_alignment2_sorted.bam &
wait
echo "Done sorting BAM file"

# Filter out mapped reads
echo "Filtering out mapped reads..."
samtools view -@ $threads2 -b -F 4 ${Base}_${Num_of_runs}_alignment1_sorted.bam > ${Base}_${Num_of_runs}_mapped_reads1.bam &
samtools view -@ $threads2 -b -F 4 ${Base}_${Num_of_runs}_alignment2_sorted.bam > ${Base}_${Num_of_runs}_mapped_reads2.bam &
wait
echo "Done filtering out mapped reads"

# Get list of mapped read names
echo "Getting list of mapped read names..."
samtools view -@ $threads2 ${Base}_${Num_of_runs}_mapped_reads1.bam | cut -f1 > ${Base}_${Num_of_runs}_mapped_read_names1.txt &
samtools view -@ $threads2 ${Base}_${Num_of_runs}_mapped_reads2.bam | cut -f1 > ${Base}_${Num_of_runs}_mapped_read_names2.txt &
wait
echo "Done getting list of mapped read names"

# Use seqtk to filter out mapped reads from the original FASTQ files
echo "Filtering out mapped reads from the original FASTQ files Read1 and Read2"
seqtk subseq $read1 ${Base}_${Num_of_runs}_mapped_read_names1.txt > ${Base}_${Num_of_runs}_filtered1.fastq &
seqtk subseq $read2 ${Base}_${Num_of_runs}_mapped_read_names2.txt > ${Base}_${Num_of_runs}_filtered2.fastq &
wait
echo "Done filtering out mapped reads from the original FASTQ files Read1 and Read2"


    # Deactivate conda
    conda deactivate || { echo "Failed to deactivate conda environment"; exit 1; }
    echo "Done moving filtered FASTA to original FASTA"
    rm ${Base}_${Num_of_runs}_mapped_read_names2.txt ${Base}_${Num_of_runs}_mapped_read_names1.txt &
    rm ${Base}_${Num_of_runs}_mapped_reads2.bam ${Base}_${Num_of_runs}_mapped_reads1.bam &
    rm ${Base}_${Num_of_runs}_alignment2_sorted.bam ${Base}_${Num_of_runs}_alignment1_sorted.bam &
    rm ${Base}_${Num_of_runs}_alignment2.bam ${Base}_${Num_of_runs}_alignment1.bam &
    rm ${Base}_${Num_of_runs}_alignment2.sam ${Base}_${Num_of_runs}_alignment1.sam &
    wait
    echo "Temp files removed"
    
    #move the re_output directory to the Repeat_explorer_outputs directory and add the number of runs to the end of the directory name
    mv ./re_output "$ProjectDir"/Repeat_explorer_outputs/re_output_${Num_of_runs} 
    # Gzip the files inside the re_output_${Num_of_runs}
    gzip "$ProjectDir"/Repeat_explorer_outputs/re_output_${Num_of_runs}/*

    # Move the filtered reads to the original read names
    mv ${Base}_${Num_of_runs}_filtered1.fastq $read1 &
    mv ${Base}_${Num_of_runs}_filtered2.fastq $read2 &
    wait
done