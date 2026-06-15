#!/bin/bash -l
#SBATCH -M calc-cpu
#SBATCH -p cpu-daily
#SBATCH -t 24:00:00
#SBATCH --job-name=MSS_MicroSensitivity
#SBATCH --mem=64G
#SBATCH --cpus-per-task=350
#SBATCH --output=logs/%x-%j.out

# Absolute path to the repository on the host
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
IMAGE_PATH="${PROJECT_DIR}/containers/mss_image.sif"
LOG_DIR="${PROJECT_DIR}/logs"

# Where to mount it inside the container
WORKDIR=/workspace

mkdir -p "${PROJECT_DIR}/logs"

singularity exec \
  --bind ${PROJECT_DIR}:${WORKDIR} \
  "${IMAGE_PATH}" \
  bash -lc "\
    cd ${WORKDIR} && \
    uv sync --frozen && \
    uv run python -m mss.cli.run_micro_sensitivity \
      --config config/cal_micro_sensitivity.yml \
  "
