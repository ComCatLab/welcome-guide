#!/bin/bash
#SBATCH --account=def-samiras
#SBATCH --job-name=Co-BHT_on_Co_vertical_1-c7zrPZuKru
#SBATCH --mem-per-cpu=1000MB
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=23:00:00
#SBATCH --mail-user=ugn@sfu.ca
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT,TIME_LIMIT_90

if test -e "/etc/profile"; then
  source "/etc/profile"
fi

if test -e "$HOME/.bash_profile"; then
  source "$HOME/.bash_profile"
fi

unset LANG

module purge
module load vasp
module load python/3.11.9

source "$COMP_CHEM_ENV"

export LC_ALL="C"
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

ulimit -s unlimited

ulimit -a -S
ulimit -a -H

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

TMP_BASE_DIR="$(dirname "$TMP_WORK_DIR")"
JOB_WORK_DIR="$(basename "$TMP_WORK_DIR")"

# Creating a symbolic link to temporary directory holding work files while job running

if ! test -e "${TMP_WORK_DIR}"; then
  mkdir "${TMP_WORK_DIR}"
fi
ln -s "${TMP_WORK_DIR}" scratch_dir
cd "${TMP_WORK_DIR}" || exit

cp -v "$SLURM_SUBMIT_DIR"/{*py,*.traj} "$TMP_WORK_DIR"/

echo " "

# Preemptively end job if getting close to time limit

export AUTOJOB_SLURM_SCRIPT="vasp.sh"

# run ase calculation and time
python3 run.py "$SLURM_SUBMIT_DIR" > "$SLURM_SUBMIT_DIR"/out.txt


AUTOJOB_FILES_TO_DELETE="EIGENVAL IBZKPT PCDAT PROCAR ELFCAR LOCPOT PROOUT TMPCAR vasp.dipcor"
rm -vf "$AUTOJOB_FILES_TO_DELETE"
sleep 10 # Sleep some time so potential stale nfs handles can disappear.

cd "${TMP_BASE_DIR}" || exit
mkdir -vp "${SLURM_SUBMIT_DIR}" # if user has deleted or moved the submit dir

tar -zcvf "${SLURM_SUBMIT_DIR}/${JOB_WORK_DIR}.tgz" "${JOB_WORK_DIR}" \
  || { echo "ERROR: Failed to create tgz-file. Please cleanup TMP_WORK_DIR $TMP_WORK_DIR on host '$HOSTNAME' manually (if not done automatically by queueing system)."; exit 102; }

rm -rvf "${TMP_WORK_DIR}"

cd "${SLURM_SUBMIT_DIR}" || exit
tar -xzf "${JOB_WORK_DIR}".tgz
mv "${JOB_WORK_DIR}"/* .
rm -r "${JOB_WORK_DIR}".tgz "${JOB_WORK_DIR}"

rm "${SLURM_SUBMIT_DIR}/scratch_dir"
