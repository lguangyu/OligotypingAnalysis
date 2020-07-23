#!/bin/bash

n_threads=4
msa=65 # consult the msa_test output
aln="mothur2oligo.fasta"
positions=$(cat $aln".oligo_pos.custom.txt")

# oligotyping needs python2
module load python/2.7.15
. /home/li.gua/python-virtual-env/python-2.7.15-metagenomic-general/bin/activate
# to generate figures it also needs R
module load R/3.6.1

entropy=$aln"-ENTROPY"
out_dir=$aln".oligotype."$positions
oligotype -M $msa -s 3 -N $n_threads -o $out_dir \
	-C $(echo $positions | sed 's/[_.]/,/g') \
	$aln $entropy

module unload R/3.6.1
deactivate
module unload python/2.7.15
