#!/bin/bash -l
#SBATCH -p cpu-daily
#SBATCH -t 24:00:00
#SBATCH --job-name=MSS_Simulation
#SBATCH --mem=64G
#SBATCH --cpus-per-task=64
#SBATCH --output=/home2/lukas.breiter/CIC_slurm_le4/MSS/logs/%x-%j.out
#SBATCH --error=/home2/lukas.breiter/CIC_slurm_le4/MSS/logs/%x-%j.err

# Absolute path to your project on the host
PROJECT_DIR=/home2/lukas.breiter/CIC_slurm_le4/MSS

# Where to mount it inside the container
WORKDIR=/workspace

mkdir -p "${PROJECT_DIR}/logs"

singularity exec \
  --bind ${PROJECT_DIR}:${WORKDIR} \
  containers/mss_image.sif \
  bash -lc "\
    cd ${WORKDIR} && \
    uv sync && \
    uv run python run_coupled_simulation.py --config shared/config_abx.yml \
  "
