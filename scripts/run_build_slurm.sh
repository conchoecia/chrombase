#!/bin/bash
#SBATCH --job-name=chrombase_build
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=10-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# Controller job for building a chrombase genome database on SLURM.
#
# Run it from the database directory (the one that holds config.yaml), with an environment
# that provides snakemake and snakemake-executor-plugin-slurm:
#
#   export CHROMBASE=/path/to/chrombase
#   export NCBI_API_KEY=...            # optional, raises the NCBI request limit
#   sbatch $CHROMBASE/scripts/run_build_slurm.sh annotated      # or: unannotated
#
# Flags that keep the database directory small and the build restartable:
#   --rerun-triggers mtime  Editing a snakefile does not re-run genomes that are already built,
#                           so there is no need to `snakemake --touch` everything after an update.
#   --drop-metadata         Do not keep a .snakemake/metadata record for every output file.
#                           Rerun decisions only compare recorded metadata, so nothing reruns
#                           because a record is missing.
#   --rerun-incomplete      Outputs of jobs that were killed are removed and rebuilt.
# The SLURM executor deletes the logs of successful jobs, so log/ only keeps failed jobs.
#
# If the controller itself is cancelled, run `snakemake --unlock` with the same snakefile
# before resubmitting.

set -euo pipefail

KIND=${1:?usage: run_build_slurm.sh annotated|unannotated}
case "$KIND" in
    annotated|unannotated) ;;
    *) echo "usage: run_build_slurm.sh annotated|unannotated" >&2; exit 1 ;;
esac

# sbatch runs a spooled copy of this script, so the repository location has to come from
# the environment when submitted that way.
if [ -z "${CHROMBASE:-}" ]; then
    if [ -n "${SLURM_JOB_ID:-}" ]; then
        echo "Set CHROMBASE to the chrombase repository before calling sbatch." >&2
        exit 1
    fi
    CHROMBASE=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
fi
SNAKEFILE="$CHROMBASE/chrombase_build_db_${KIND}_chr.snakefile"
[ -f "$SNAKEFILE" ] || { echo "Snakefile not found: $SNAKEFILE" >&2; exit 1; }
[ -f config.yaml ] || { echo "Run this from the database directory (no config.yaml in $PWD)." >&2; exit 1; }

if [ -z "${NCBI_API_KEY:-}" ]; then
    echo "NCBI_API_KEY is not set; NCBI requests will use the lower anonymous rate limit." >&2
fi

snakemake --snakefile "$SNAKEFILE" \
    --executor slurm \
    --jobs "${JOBS:-700}" \
    --resources download_slots="${DOWNLOAD_SLOTS:-25}" \
    --default-resources mem_mb=8000 runtime=60 \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --drop-metadata \
    --slurm-logdir log
