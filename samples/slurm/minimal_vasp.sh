#!/bin/bash
#SBATCH --account=def-samiras-ab
#SBATCH --job-name=JOB_NAME
#SBATCH --mem-per-cpu=1000MB
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=23:00:00
#SBATCH --mail-user=SFU_ID@sfu.ca
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT,TIME_LIMIT_90

# load personal configuration
if test -e "/etc/profile"; then source "/etc/profile"; fi
if test -e "$HOME/.bash_profile"; then source "$HOME/.bash_profile"; fi

# software setup
unset LANG
ulimit -s unlimited
setopt enable aliases
export LC_ALL="C"
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
module purge
module load vasp
activate_env

# run ase calculation
python3 "$AUTOJOB_PYTHON_SCRIPT"

# record job in log file
echo "${SLURM_JOB_ID}-${SLURM_JOB_NAME}" is complete: on "$(date +'%y.%m.%d %H:%M:%S')" "${SLURM_SUBMIT_DIR}" >> ~/job.log
