#!/bin/bash
#SBATCH --job-name=Ectopic_recombination
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jlamb@crimson.ua.edu
#SBATCH -N 1 #nodes
#SBATCH -o Ectopic_%j.o

#set up
echo "Starting the pipeline at $(date)"
echo "Setting up the environment"
threads=32
#divide the number of threads by 2
threads2=$((threads / 2))
#set the species hardcode for now
#in the future it will be passed as an argument
Base="SRX19958881"
echo "Processing samples with base name: $Base"

ProjectDir="/home/jlamb1/Projects/Ectopic/$Base"
[ ! -d "$ProjectDir" ] && mkdir -p "$ProjectDir"
cd "$ProjectDir" || { echo "Failed to change directory to $ProjectDir"; exit 1; }
mkdir -p "$ProjectDir"/Repeat_explorer_outputs

echo "Copying files to Directory"
start=$(date +%s.%N)
cp /home/jlamb1/SRA_seqs/"$Base"_N1.fastq.gz /home/jlamb1/Projects/Ectopic/"$Base"/ &
cp /home/jlamb1/SRA_seqs/"$Base"_N2.fastq.gz /home/jlamb1/Projects/Ectopic/"$Base"/ &
wait
end=$(date +%s.%N)
runtime=$(echo "$end - $start" | bc)
echo "Done copying files to Directory in $runtime seconds"

#if /home/jlamb1/SRA_seqs/"$Base"_N1 and N2.fastq exists, then remove the file
if [ -f /home/jlamb1/Projects/Ectopic/"$Base"/"$Base"_N1.fastq ]; then
rm /home/jlamb1/Projects/Ectopic/"$Base"/"$Base"_N1.fastq
fi
if [ -f /home/jlamb1/Projects/Ectopic/"$Base"/"$Base"_N2.fastq ]; then
rm /home/jlamb1/Projects/Ectopic/"$Base"/"$Base"_N2.fastq
echo "Removed existing fastq files"
fi

#gunzip the files if they end in .gz and overwrite the files
echo "Unzipping files"
start=$(date +%s.%N)
gunzip -f /home/jlamb1/Projects/Ectopic/"$Base"/"$Base"_N1.fastq.gz &
gunzip -f /home/jlamb1/Projects/Ectopic/"$Base"/"$Base"_N2.fastq.gz &
wait
end=$(date +%s.%N)
runtime=$(echo "$end - $start" | bc)
echo "Done unzipping files to Directory in $runtime seconds"

read1="/home/jlamb1/Projects/Ectopic/$Base/${Base}_N1.fastq"
read2="/home/jlamb1/Projects/Ectopic/$Base/${Base}_N2.fastq"

echo "Finished setting up the environment at $(date)"


