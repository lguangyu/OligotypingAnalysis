#!/bin/bash

n_threads=4
msa=65
n_pos=7
aln="mothur2oligo.fasta"

# oligotyping needs python2
module load python/2.7.15
. /home/li.gua/python-virtual-env/python-2.7.15-metagenomic-general/bin/activate
# to generate figures it also needs R
module load R/3.6.1

entropy=$aln"-ENTROPY"
out_dir=$aln".position_oligotype"
# find the positions
positions=$(sort -nrk2 $entropy | head -n $n_pos | cut -f1 |\
	perl -ne 'chomp; print($_, ",");' | sed 's/,$//')
echo $positions
oligotype -M $msa -s 3 -C $positions -N $n_threads -o $out_dir \
	$aln $aln"-ENTROPY"

module unload R/3.6.1
deactivate
module unload python/2.7.15
