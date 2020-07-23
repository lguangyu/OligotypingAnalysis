#!/bin/bash

input_fasta="mothur2oligo.fasta"

# oligotyping needs python2
module load python/2.7.15
. /home/li.gua/python-virtual-env/python-2.7.15-metagenomic-general/bin/activate
entropy-analysis --no-display $input_fasta
deactivate
module unload python/2.7.15
