#!/usr/bin/env bash
snakemake -k --printshellcmds -j 500 --cluster "sbatch -J OrthoMam{params.name} -p {params.queue} -N 1 -o ./slurm/%x.%j.out -e ./slurm/%x.%j.err --cpus-per-task={params.threads} --mem={params.mem} -t {params.time}"
