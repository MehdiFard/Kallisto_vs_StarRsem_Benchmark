#!/bin/bash
### kallisto vs STAR+RSEM benchmark - Snakemake controller job (submits every rule as its own SLURM job).
### Edit the #SBATCH account/partition lines and profiles/slurm/config.yaml, then submit from the project root:
###   sbatch run_pipeline.slurm                                   full pipeline (target: all)
###   sbatch --export=ALL,TARGET=quant run_pipeline.slurm         quantification + comparison tables only
###   sbatch --export=ALL,TARGET=trans run_pipeline.slurm         optional trans-eQTL
###   sbatch --export=ALL,UNLOCK=1 run_pipeline.slurm             after a killed controller (removes its lock)
###   sbatch --export=ALL,SNAKEMAKE_ARGS="--config samples=config/samples_pilot.tsv",TARGET=quant run_pipeline.slurm
### SNAKEMAKE_ENV: conda env with snakemake >= 8 and snakemake-executor-plugin-slurm (default: snakemake).

#SBATCH --job-name=kallisto_vs_star
#SBATCH --output=logs/slurm/%x-%j.log
#SBATCH --error=logs/slurm/%x-%j.err
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_CPU_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00

eval "$(conda shell.bash hook)"
conda activate "${SNAKEMAKE_ENV:-snakemake}"
set -euo pipefail   # after activation: some conda activate.d scripts use unset variables

PROFILE=${PROFILE:-profiles/slurm}
if [ "${UNLOCK:-0}" = 1 ]; then snakemake --profile "$PROFILE" --unlock; fi
snakemake --profile "$PROFILE" ${TARGET:-all} ${SNAKEMAKE_ARGS:-}
