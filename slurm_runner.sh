#!/bin/bash -l
#SBATCH -p cpu-daily
#SBATCH -t 24:00:00
#SBATCH --job-name=MSS_Simulation
#SBATCH --mem=64G
#SBATCH --cpus-per-task=64
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Absolute path to the repository on the host
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE_PATH="${PROJECT_DIR}/containers/mss_image.sif"

# Where to mount it inside the container
WORKDIR=/workspace

mkdir -p "${PROJECT_DIR}/logs"

if [ ! -f "${PROJECT_DIR}/pyproject.toml" ]; then
  echo "Expected pyproject.toml in ${PROJECT_DIR}, but it was not found." >&2
  exit 1
fi

if [ ! -f "${IMAGE_PATH}" ]; then
  echo "Expected Singularity image at ${IMAGE_PATH}, but it was not found." >&2
  exit 1
fi

singularity exec \
  --bind ${PROJECT_DIR}:${WORKDIR} \
  "${IMAGE_PATH}" \
  bash -lc "\
    cd ${WORKDIR} && \
    uv sync && \
    uv run python run_coupled_simulation.py --config shared/config_abx.yml \
  "
