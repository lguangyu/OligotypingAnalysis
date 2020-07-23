#!/bin/bash

ln -sfT ../mothur ./mothur
ln -sfT ./mothur/mothur.output ./mothur.output
ln -sfT ./mothur.output/mothur.input.trim.contigs.good.unique.good.filter.denovo.vsearch.pick.pick.count_table ./mothur.output.count_table
ln -sfT ./mothur.output/mothur.input.trim.contigs.good.unique.good.filter.unique.pick.pick.fasta ./mothur.output.fasta
ln -sfT ./mothur.output/mothur.input.trim.contigs.good.unique.good.filter.unique.pick.classif.wang.pick.taxonomy ./mothur.output.tax
