# Wrapped Oligotyping Analysis Pipeline

## Dependencies

* Unix-like environment with QoL tools like `grep` and `sed`
* mothur 1.48.0 (required by mothur SOP)
* python3
* oligotyping (python3 package) and its dependencies
* perl (as needed by `mothur2oligo` scripts)
* mafft (to realign sequences)
* (optional) if you want oligotyping's native plot features, R and ggplot2 are required

## Installation guide

### Install python package dependencies

All dependency python packages can be installed by

```bash
$ pip install -r requirements.txt
```

### Install oligotyping

The original package of `oligotyping` contains an error to halt automated `pip` installation. We need to resolve the error manually. First download the source package from `pip`:

```bash
$ mkdir temp && cd temp
$ pip3 download --no-deps oligotyping
$ tar -zxf oligotyping-3.1.tar.gz
$ cd oligotyping-3.1
```

This will download only `oligotyping` without all its python dependencies, and extract contents from the tarball. Then edit the file `bin/o-boxplots.R` at lines 75, 84, 97, and 125, delete the weird characters as shown below:

```
    #Êperform lda on the node with respect to classes
     ^
```

Last, install with `pip`:

```bash
$ pip install ./
```

## Before oligotyping analysis

Run the mothur SOP in `mothur/`, following the procedure described in `mothur/README.md`.

## Oligotyping example run

### 1. Copy analysis template

Oligotyping can only be done at per-taxon basis, as it is meant to resolve the OTUs of an interested taxon beyond genus level.

The `oligo.prototype` is a analysis template for each targeted taxon. The analysis begins with cloning the `oligo.prototype` directory, e.g. one copy for each taxon targetted for analysis. Here we'll use the genus `Acinetobacter` as an exampe:

```bash
$ copy -r oligo.prototype/ oligo.acinetobacter/
```

### 2. Find the full taxonomy

It is essential to correctly find the full name of the taxonomy for analysis. This information can be acquired from the `mothur` output. We can do a quick `grep` search:

```bash
# get into the mothur2oligo directory
$ cd oligo.acinetobacter/mothur2oligo
$ grep -i 'acinetobacter' mothur.output.seqs.taxonomy
```

Note the file `mothur.output.seqs.taxonomy` is a symbolic link to the actual mothur output file. If you setup the environment correctly and `mothur` SOP finished correctly, these files should exist.

Output of the above `grep` command may look like something below (showing first 5 lines):

```
ASV0001	4643	Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100);Acinetobacter_unclassified(100);
ASV0007	1160	Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100);midas_s_32693(100);
ASV0012	569	Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100);Acinetobacter_parvus(100);
ASV0104	49	Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100);Acinetobacter_unclassified(100);
ASV0166	28	Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100);Acinetobacter_towneri(100);
```

It is a 3-column table. What we need is in the 3rd column, as the taxonomy path from root to the genus level of `Acinetobacter`, which is:

```
Bacteria(100);Proteobacteria(100);Gammaproteobacteria(100);Pseudomonadales(100);Moraxellaceae(100);Acinetobacter(100);
```

Then discard all the bootstrapping numbers (with the parenthese). The final taxonomy string should look like something below:

```
Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Moraxellaceae;Acinetobacter;
```

Note that we need to preserve every single semicolun `;` including the last one. Then save the above taxonomy text in the file `extract_taxon` (using any text editor).

### 3. Adapt `mothur` output for `oligotyping` analysis.

This is done using the `script/mothur2oligo.sh` script. The procedure contains:

1. extract all unique sequences belong to the taxonomy in `extract_taxon`
2. realign using `mafft`
3. dereplicate the sequences
4. rename the sequences to make `oligotyping` happy

To do all above, run:

```bash
$ bash script/mothur2oligo.sh
```

Note that this script is designed to work with `slurm` in a cluster-like environment. You might need to configure it accordingly (e.g. configure correct paths to `mothur`/`mafft` and number of processors) to run it on your own machine or cluster.

### 4. Entropy analysis

Now we are done with `mothur2oligo`, and need to swtich to `oligotyping`:

```bash
$ cd ../oligotyping
```

Then run:

```bash
$ bash script/entropy_analysis.sh
```

It will generate abunch of files, the most important are:

* `mothur2oligo.fasta-ENTROPY`: the file with entropy values at each location
* `filtered_positions`: the positions list with an entropy values >= 0.2

### 5. Oligotyping

Now we run oligotyping based on the entropy just calculated and the positions selected:

```bash
$ bash script/oligotyping.sh
```

The results will be saved in `mothur2oligo.fasta.oligo_final`.

### 6. Oligotype filtering

Note the above approach is a very loose approach that will result in a lot of oligotypes apparently. This is because we haven't done any filtering yet. The "official" way to do is to pass the `-M` argument a positive integer to the `oligotype` script in `script/oligotyping.sh`. This argument will filter out any oligos that have a count number lower it. However I find this hard-coded way is not flexible and requires some human intervention as the threashold may change based on both the abundance of the targeted taxon and sequncing depth. Alternatively, I decide to go another way using custom script, for example:

```bash
$ python script/get_abund_oligo_list.py \
	mothur2oligo.fasta.oligo_final \
	--abund-threshold 0.1 \
	--count-threshold 10 \
	-o abund_oligo.list
```

This script filters oligos if:

* is more abundant (w.r.t. extracted taxon=100%) than 0.1 (10%) in at least one sample, AND:
* has more than 10 counts in total (this functionality reprecates the previous `-M`)

Note that you can set either argument to 0 (or omit the argument) to disable the corresponding filter, for example:

```bash
--abund-threshold 0.1
--count-threshold 0
```

will only do abundance threshold filtering but not the count threshold filtering; This is exactly the same as running:

```bash
--abund-threshold 0.1
```

### 7. Downstream analysis

Up to this point, the oligotyping analysis is done. The essential output files are:

* `mothur2oligo.final.oligo_final/MATRIX-COUNT.txt`: the table of oligo count numbers in each sample
* `mothur2oligo.final.oligo_final/MATRIX-PERCENT.txt`: the table of oligo abundance percentages in each sample; this is essentially the normalized version of `MATRIX-COUNT.txt`
* `abund_oligo.list`: the list of filtered oligos

There are more things can be interesting, for example determining the taxonomy of each oligo. Those are considered downstream analysis. Since the approaches are many, they will not be included in this example. One possible approach is to exhausively search the taxonomy classification of every sequences in an oligotype (do not use the representative sequences) against NCBI's RNA refseq database then determine the oligotype taxonomy via majority vote. However considering the number of oligotypes and the size of database, it must be done with HPC.