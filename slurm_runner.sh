#!/bin/bash -l
#SBATCH -p cpu-daily
#SBATCH -t 24:00:00
#SBATCH --job-name=MSS_Simulation
#SBATCH --mem=64G
#SBATCH --cpus-per-task=64

# Absolute path to your project on the host
PROJECT_DIR=/home2/lukas.breiter/CIC_slurm_le4/MSS

# Where to mount it inside the container
WORKDIR=/workspace

singularity exec \
     --bind ${PROJECT_DIR}:${WORKDIR} \
     python_poetry.sif \
     bash -lc "\
       cd ${WORKDIR} && \
       uv sync --no-interaction && \
       uv run python run_coupled_simulation.py --config shared/config_abx.yml \
     "
