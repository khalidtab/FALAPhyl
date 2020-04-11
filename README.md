# What does this snakemake workflow do?
Automatic exploratory analysis of microbiome data.

# How do you run it?
For the workflow to run, snakemake expects the names of the output files it is supposed to generate. This workflow has been made with the idea that each file name (without its original extension) is to be passed to snakemake with the extension of ".final". For example, "smokers_nonsmokers.final". Of course you don't have to pass these values interactively to snakemake and the values can be saved to a text file then fed to snakemake

> samples=$(cat samples.txt) &&  snakemake $samples --cores all --use-conda

# Necessary folder structure

Before execution, the snakemake folder needs the following folders inside of it where the biom files and mapping files reside in:

- data/biom
- data/map 

Each biom and mapping file needs to be named the same, and be stored in the correct folders. For example: 

- data/biom/health_disease.biom
- data/map/health_disease.txt