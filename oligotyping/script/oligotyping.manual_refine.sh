#!/bin/bash

n_threads=4
msa=65
man_pos1="18,81,95,102,268,476,491,508"
positions="83,245,105,155,185,107,244,"$man_pos1
aln="mothur2oligo.fasta"

# oligotyping needs python2
module load python/2.7.15
. /home/li.gua/python-virtual-env/python-2.7.15-metagenomic-general/bin/activate
# to generate figures it also needs R
module load R/3.6.1

entropy=$aln"-ENTROPY"
out_dir=$aln".position_oligotype."$(echo $man_pos1 | sed 's/,/_/g')
# find the positions
oligotype -M $msa -s 3 -C $positions -N $n_threads -o $out_dir \
	$aln $aln"-ENTROPY"

module unload R/3.6.1
deactivate
module unload python/2.7.15
