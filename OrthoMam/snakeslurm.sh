#!/usr/bin/env bash
snakemake --forcerun -k --printshellcmds -j 512 --cluster "sbatch -J OrthoMam{params.name} -p cpu -N 1 -o ./slurm/%x.%j.out -e ./slurm/%x.%j.err --cpus-per-task={params.threads} --mem={params.mem} -t {params.time}"
