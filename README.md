## Overview
The mustache workflow is designed to identify insertion events from paired-end metagenomic sequencing data. The required inputs are FASTQ files, reference genomes for bacterial species of interest, reference sequences for insertion sequences of interst, and read-level taxonomic classifications for both forward and reverse reads. The approach if this workflow enables scientists to identify longitudinal changes in strain-level insertion event frequency, as well as variations between individuals and body sites.

## Installation

mustache is implemented as a snakemake workflow. Snakemake is a workflow management system that allows for all the tasks to be easily parallelized. I strongly recommend you read about snakemake and even complete the beginner's tutorial [here](https://snakemake.readthedocs.io/en/stable/).
