snakemake -s setup_run.smk --profile slurm
snakemake -s run_extractor.smk --profile slurm
snakemake -s finish_process.smk --profile slurm
