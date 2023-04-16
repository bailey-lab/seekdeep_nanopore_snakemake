# seekdeep_nanopore
a basic workflow for running seekdeep on nanopore

## Installation:
Install conda with:
https://github.com/conda-forge/miniforge#mambaforge

Install snakemake in an environment called snakemake with:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```

## Usage:
 - Change directory to a folder where you want to run the analysis
 - Download the seekdeep_nanopore_general.smk file into this folder
 - Download the sif file from here into the same folder: https://seekdeep.brown.edu/programs/elucidator.sif
 - Download the seekdeep_nanopore_general.yaml file into the same folder
 - Edit the config.yaml file using the instructions in the comments. Use a text editor that outputs unix line endings (e.g. vscode, notepad++, gedit, micro, emacs, vim, vi, etc.)
 - Activate snakemake with:
```bash
conda activate snakemake
```
 - Run snakemake with:
```bash
snakemake -s seekdeep_nanopore_general.smk --cores [your_desired_core_count]
```
