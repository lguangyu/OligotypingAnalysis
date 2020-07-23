#!/bin/bash

aln="mothur2oligo.fasta"
entropy=$aln"-ENTROPY"
output=$aln".oligo_pos.custom.txt"
n_pos=20 # check entropy output

# print
sort -nrk2 $entropy |\
	head -n $n_pos
# make output
sort -nrk2 $entropy |\
	head -n $n_pos |\
	cut -f 1 |\
	perl -ne 'chomp; print($_, "_");' |\
	sed 's/_$//' |\
	tee $output
