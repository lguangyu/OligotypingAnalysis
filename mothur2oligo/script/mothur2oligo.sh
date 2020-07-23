# USAGE: sh mothur2oligo.sh
# This is a shell script for transforming mothur output to appropriate format for 
# A. murat Eren's oligotyping pipeline 

## Set variables

# Adjust the file names to your own study - these are the files from the mothur SOP
taxonomy="mothur.output.tax"
fasta="mothur.output.fasta"
count="mothur.output.count_table"
processors=1

# Set the taxon you want to select, separate taxonomic levels with ";" 
# Do not touch inner and outer quotes
taxon="'Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Moraxellaceae;Acinetobacter;'"


################################
########## Script  #############
################################

redundantFasta=$(echo ${fasta}.pick.redundant.fasta | sed 's/.fasta//')
groups=$(echo ${count}.pick.redundant.groups | sed 's/.count_table//') 

# Call mothur commands for generating deuniqued fasta file for a specific lineage
mothur "#set.current(processors=$processors); get.lineage(taxonomy=$taxonomy, taxon=$taxon, count=$count); list.seqs(count=current); get.seqs(accnos=current, fasta=$fasta); deunique.seqs(fasta=current, count=current)"

# Replace all "_" in fasta header with a ":"
cat $groups | sed 's/_/:/g' > intermediate1
# Make a file which maps sample names to sequence headers
paste $groups intermediate1 | awk 'BEGIN{FS="\t"}{print $1"\t"$2"_"$3}' > intermediate2

# Perl script to rename the headers of the fasta to include the sample name at the beginning followed by a "_"
perl ./script/renamer.pl $redundantFasta intermediate2 

