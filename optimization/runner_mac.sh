#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OPTIMIZATION_SCRIPT="${SCRIPT_DIR}/optimization.py"
CONDA_ENV_NAME="${UTAG_CONDA_ENV:-microlive}"

usage() {
    cat <<'EOF'
Usage:
  ./runner_mac.sh [dataset] [workers]

Datasets:
  0 or utag          UTag (default)
  1 or utag_c_free   UTag C-free
  2 or suntag        SunTag
  3 or alfatag       ALFATag

Examples:
  ./runner_mac.sh utag
  ./runner_mac.sh suntag 6
  ./runner_mac.sh 3 8

Optional environment variable:
  UTAG_CONDA_ENV      Conda environment name (default: microlive)
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

case "${1:-0}" in
    0|utag) dataset_index=0; dataset_name="utag" ;;
    1|utag_c_free|utag-c-free) dataset_index=1; dataset_name="utag_c_free" ;;
    2|suntag) dataset_index=2; dataset_name="suntag" ;;
    3|alfatag) dataset_index=3; dataset_name="alfatag" ;;
    *)
        echo "Error: unknown dataset '${1}'." >&2
        usage >&2
        exit 2
        ;;
esac

logical_cpus="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 4)"
default_workers=$((logical_cpus > 2 ? logical_cpus - 2 : 1))
if ((default_workers > 8)); then
    default_workers=8
fi
workers="${2:-${default_workers}}"

if ! [[ "${workers}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: workers must be a positive integer." >&2
    exit 2
fi

if ((workers > logical_cpus)); then
    echo "Error: requested ${workers} workers, but this Mac reports ${logical_cpus} logical CPUs." >&2
    exit 2
fi

if [[ ! -f "${OPTIMIZATION_SCRIPT}" ]]; then
    echo "Error: optimization.py not found at ${OPTIMIZATION_SCRIPT}" >&2
    exit 1
fi

acf_file="${SCRIPT_DIR}/source_files_ACF/data_ACF/df_ACF_${dataset_name}.csv"
if [[ ! -f "${acf_file}" ]]; then
    echo "Error: input ACF file not found at ${acf_file}" >&2
    exit 1
fi

if command -v conda >/dev/null 2>&1; then
    conda_executable="$(command -v conda)"
elif [[ -x /opt/anaconda3/bin/conda ]]; then
    conda_executable="/opt/anaconda3/bin/conda"
else
    echo "Error: Conda was not found. Install Anaconda/Miniconda or add conda to PATH." >&2
    exit 1
fi

conda_base="$("${conda_executable}" info --base)"
source "${conda_base}/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -Fxq "${CONDA_ENV_NAME}"; then
    echo "Error: Conda environment '${CONDA_ENV_NAME}' does not exist." >&2
    echo "Create it with: conda env create -f ${PROJECT_ROOT}/environment.yml" >&2
    exit 1
fi

conda activate "${CONDA_ENV_NAME}"

task_tmp_root="${TMPDIR:-/tmp}/utag_optimization"
export NUMBA_CACHE_DIR="${task_tmp_root}/numba_cache"
export MPLCONFIGDIR="${task_tmp_root}/matplotlib_cache"
export UTAG_N_JOBS="${workers}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "${NUMBA_CACHE_DIR}" "${MPLCONFIGDIR}"

results_dir="${SCRIPT_DIR}/results_ACF"
logs_dir="${results_dir}/logs"
mkdir -p "${logs_dir}"
timestamp="$(date '+%Y%m%d_%H%M%S')"
log_file="${logs_dir}/optimization_${dataset_name}_${timestamp}.log"

echo "Project:     ${PROJECT_ROOT}"
echo "Dataset:     ${dataset_name} (index ${dataset_index})"
echo "Environment: ${CONDA_ENV_NAME}"
echo "Workers:     ${workers} of ${logical_cpus} logical CPUs"
echo "Inputs:      ${SCRIPT_DIR}/source_files_ACF"
echo "Results:     ${results_dir}"
echo "Log:         ${log_file}"
echo
echo "The optimizer will reuse compatible files in simulation_cache."
echo "Press Control-C to stop. Completed cache files can be reused on the next run."
echo

cd "${PROJECT_ROOT}"

if command -v caffeinate >/dev/null 2>&1; then
    caffeinate -i python -u "${OPTIMIZATION_SCRIPT}" "${dataset_index}" 2>&1 | tee "${log_file}"
else
    python -u "${OPTIMIZATION_SCRIPT}" "${dataset_index}" 2>&1 | tee "${log_file}"
fi
