#!/bin/bash

n_threads=4
positions="83,245,105,155,185,107,244"
aln="mothur2oligo.fasta"

# oligotyping needs python2
module load python/2.7.15
. /home/li.gua/python-virtual-env/python-2.7.15-metagenomic-general/bin/activate

# remove old results
msa_dir=$aln".msa_test"
rm -rf $msa_dir
mkdir $msa_dir
res_file=$aln".msa_test.tsv"
rm -f $res_file

# the parameter value to test
msa_list=$msa_dir"/msa.list"
seq 30 2 120 > $msa_list

# iterative test -M
entropy=$aln"-ENTROPY"
for msa in $(cat $msa_list); do
    # output directory for each run
    out_dir=$msa_dir"/"$msa
    # run test
    oligotype -M $msa -s 3 -C $positions -N $n_threads -o $out_dir \
        $aln $entropy
    # extract results
    res=$(grep -Ee '-M.*elimination' $out_dir"/HTML-OUTPUT/index.html" |\
        grep -oP '\d+\D*$' | grep -oP '\d+')
    echo -e "$msa\t$res" >> $res_file
done

deactivate
module unload python/2.7.15
