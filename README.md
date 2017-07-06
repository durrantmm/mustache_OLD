# mustache
The mustache workflow is designed to identify insertion events from paired-end metagenomic sequencing data. The required inputs are FASTQ files, reference genomes for bacterial species of interest, reference sequences for insertion sequences of interst, and read-level taxonomic classifications for both forward and reverse reads. The approach if this workflow enables scientists to identify longitudinal changes in strain-level insertion event frequency, as well as variations between individuals and body sites.

## Installation
Mustache is implemented as a snakemake workflow. Snakemake is a workflow management system that allows for all the tasks to be easily parallelized. I strongly recommend you read about snakemake and even complete the beginner's tutorial [here](https://snakemake.readthedocs.io/en/stable/).

This pipeline requires you to create a conda environment, so that you can then easily download all of the required packages and run the workflow on any system you'd like. If you haven't already, download anaconda [here](https://www.continuum.io/downloads). I strongly suggest you learn the basics of anaconda before continuing with installation.

### Download mustache
From a unix terminal, type the following:

~~~~
git clone https://github.com/durrantmm/mustache.git
cd mustache
~~~~

 You are now in the downloaded `mustache` directory.
 
 Let's now create a conda environment from the provided `environment.yaml` file in the `mustache` directory:

~~~~
conda env create -f environment.yaml
~~~~

 When asked if you want to install of the packages, type `y` and enter.


