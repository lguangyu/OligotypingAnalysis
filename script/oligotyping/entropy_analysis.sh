#!/bin/bash
#SBATCH -J OLIGO_ENTROPY_ANALYSIS
#SBATCH -pshort -N1 -c1

input_fasta="mothur2oligo.fasta"

entropy-analysis --no-display $input_fasta

# post process
sort -rnk2 ${input_fasta}-ENTROPY > ${input_fasta}-ENTROPY.ranked

script/filter_position.py \
	-t 0.2 \
	-o filtered_positions \
	${input_fasta}-ENTROPY.ranked
