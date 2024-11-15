#!/bin/bash
#SBATCH --account=def-samiras
#SBATCH --job-name=JOB_NAME
#SBATCH --mem-per-cpu=1000MB
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=23:00:00
#SBATCH --mail-user=SFU_ID@sfu.ca
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT,TIME_LIMIT_90

echo " "
echo "### Setting up shell environment ..."
echo " "

if test -e "/etc/profile"; then
  source "/etc/profile"
fi

if test -e "$HOME/.bash_profile"; then
  source "$HOME/.bash_profile"
fi

unset LANG

module purge
module load gaussian/g16.c01
module load python/3.11.9

# Replace "$COMP_CHEM_ENV" with the path to your Python virtual environment
source "$COMP_CHEM_ENV"

export LC_ALL="C"
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

ulimit -s unlimited

echo " "
echo "### Printing basic job infos to stdout ..."
echo " "
echo "START_TIME               = $(date '+%y-%m-%d %H:%M:%S %s')"
echo "HOSTNAME                 = ${HOSTNAME}"
echo "USER                     = ${USER}"
echo "SLURM_JOB_NAME           = ${SLURM_JOB_NAME}"
echo "SLURM_JOB_ID             = ${SLURM_JOB_ID}"
echo "SLURM_SUBMIT_DIR         = ${SLURM_SUBMIT_DIR}"
echo "SLURM_JOB_NUM_NODES      = ${SLURM_JOB_NUM_NODES}"
echo "SLURM_NTASKS             = ${SLURM_NTASKS}"
echo "SLURM_NODELIST           = ${SLURM_NODELIST}"
echo "SLURM_JOB_NODELIST       = ${SLURM_JOB_NODELIST}"

if test -f "${SLURM_JOB_NODELIST}"; then
  echo "SLURM_JOB_NODELIST (begin) ----------"
  cat "${SLURM_JOB_NODELIST}"
  echo "SLURM_JOB_NODELIST (end) ------------"
fi

echo "--------------- ulimit -a -S ---------------"
ulimit -a -S
echo "--------------- ulimit -a -H ---------------"
ulimit -a -H
echo "----------------------------------------------"

echo " "
echo "### Creating TMP_WORK_DIR directory and changing to it ..."
echo " "

if test -e "$HOME/scratch"; then
  TMP_WORK_DIR="$HOME/scratch/${SLURM_JOB_ID}"
elif test -e /scratch/"${SLURM_JOB_ID}"; then
  TMP_WORK_DIR=/scratch/${SLURM_JOB_ID}
else
  TMP_WORK_DIR="$(pwd)"
fi

# Pass memory request, cpu list, and scratch directory to Gaussian
export GAUSS_MDEF="${SLURM_MEM_PER_NODE}MB"
GAUSS_CDEF=$(taskset -cp $$ | awk -F':' '{print $2}')
export GAUSS_CDEF
export GAUSS_SCRDIR=${TMP_WORK_DIR}

TMP_BASE_DIR="$(dirname "$TMP_WORK_DIR")"
JOB_WORK_DIR="$(basename "$TMP_WORK_DIR")"

echo "TMP_WORK_DIR         = ${TMP_WORK_DIR}"
echo "TMP_BASE_DIR         = ${TMP_BASE_DIR}"
echo "JOB_WORK_DIR         = ${JOB_WORK_DIR}"

# Creating a symbolic link to temporary directory holding work files while job running

ln -s "${TMP_WORK_DIR}" scratch_dir
cd "${TMP_WORK_DIR}" || exit

echo " "
echo "### Copying input files for job (if required):"
echo " "

script_name="${BASH_SOURCE[0]}"
AUTOJOB_SLURM_SCRIPT="$(basename "$script_name")"
export AUTOJOB_SLURM_SCRIPT
export AUTOJOB_PYTHON_SCRIPT="run.py"
export AUTOJOB_COPY_TO_SCRATCH="*.chk,*.py,*.traj,*.rwf"
cp -v "$SLURM_SUBMIT_DIR"/{*.chk,*.py,*.traj,*.rwf} "$TMP_WORK_DIR"/

echo " "

# Preemptively end job if getting close to time limit
timeline=$(grep -E -m 1 '^#SBATCH[[:space:]]*--time=' "$script_name")
timeslurm=${timeline##*=}
IFS=- read -ra day_split_time <<< "$timeslurm"
no_days_time=${day_split_time[1]}
days=${no_days_time:+${day_split_time[0]}}
no_days_time=${day_split_time[1]:-${day_split_time[0]}}
IFS=: read -ra split_time <<< "$no_days_time"

# Time formats with days: D-H, D-H:M, D-H:M:S
if [[ $days ]]; then
  slurm_days="$days"
  slurm_hours=${split_time[0]}
  slurm_minutes=${split_time[1]:-0}
  slurm_seconds=${split_time[2]:-0}
# Time format without days: M, M:S, H:M:S
else
  slurm_days=0
  if [[ ${#split_time[*]} == 3 ]]; then
    slurm_hours=${split_time[0]}
    slurm_minutes=${split_time[1]}
    slurm_seconds=${split_time[2]}
  else
    slurm_hours=0
    slurm_minutes=${split_time[0]}
    slurm_seconds=${split_time[1]:-0}
  fi
fi

echo "Running for $(echo "$slurm_days*1" |bc)d $(echo "$slurm_hours*1" |bc)h $(echo "$slurm_minutes*1" |bc)m and $(echo "$slurm_seconds*1" |bc)s."
timeslurm=$(echo "$slurm_days*86400 + $slurm_hours*3600 + $slurm_minutes*60 + $slurm_seconds" | bc)
echo "This means $timeslurm seconds."
timeslurm=$(echo "$timeslurm *0.9" |bc)
echo "Will terminate at ${timeslurm}s to copy back necessary files from scratch"

echo ""
echo ""

# run ase calculation and time
timeout "${timeslurm}" python3 "$AUTOJOB_PYTHON_SCRIPT"

exit_code=$?

if [ "$exit_code" -eq 124 ]; then
  echo " "
  echo "Cancelled due to time limit."
else
  echo " "
  echo "Time limit not reached."
fi

echo " "
echo "### Cleaning up files ... removing unnecessary scratch files ..."
echo " "

AUTOJOB_FILES_TO_DELETE="*.d2e *.int *.rwf *.skr *.inp"
rm -vf "$AUTOJOB_FILES_TO_DELETE"
sleep 10 # Sleep some time so potential stale nfs handles can disappear.

echo " "
echo "### Compressing results and copying back result archive ..."
echo " "
cd "${TMP_BASE_DIR}" || exit
mkdir -vp "${SLURM_SUBMIT_DIR}" # if user has deleted or moved the submit dir

echo " "
echo "Creating result tgz-file '${SLURM_SUBMIT_DIR}/${JOB_WORK_DIR}.tgz' ..."
echo " "

tar -zcvf "${SLURM_SUBMIT_DIR}/${JOB_WORK_DIR}.tgz" "${JOB_WORK_DIR}" \
  || { echo "ERROR: Failed to create tgz-file. Please cleanup TMP_WORK_DIR $TMP_WORK_DIR on host '$HOSTNAME' manually (if not done automatically by queueing system)."; exit 102; }

echo " "
echo "### Remove TMP_WORK_DIR ..."
echo " "
rm -rvf "${TMP_WORK_DIR}"

echo " "
echo "Extracting result tgz-file"
echo " "

cd "${SLURM_SUBMIT_DIR}" || exit
tar -xzf "${JOB_WORK_DIR}".tgz
mv "${JOB_WORK_DIR}"/* .
rm -r "${JOB_WORK_DIR}".tgz "${JOB_WORK_DIR}"

rm "${SLURM_SUBMIT_DIR}/scratch_dir"
echo "END_TIME             = $(date +'%y-%m-%d %H:%M:%S %s')"

# Record job in log file

echo "${SLURM_JOB_ID}-${SLURM_JOB_NAME}" is complete: on "$(date +'%y.%m.%d %H:%M:%S')" "${SLURM_SUBMIT_DIR}" >> ~/job.log

echo " "
echo "### Exiting with exit code ${exit_code}..."
echo " "
exit "$exit_code"
