#!/bin/bash
# This script simply initiates an interactive session on the cluster to avoid using the login node.
# The session has 4 cores, 12 GB memory, X11 forwarding, and a 1 hour time limit.
# SBATCH --x11
#SBATCH --time=00:10:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=3G
#SBATCH --account=$SBATCH_ACCOUNT
#SBATCH --mail-user=[Alliance-linked-email-address
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT,TIME_LIMIT_90,FAIL

# This loads the preliminary modules necessary for VASP, then modules necessary for ASE environment, and finally loads VASP
module purge
module load StdEnv/2023 intel/2023.2.1 openmpi/4.1.5
module load python/3.11.5 scipy-stack
module load vasp/6.4.2

# This checks if the modules are loaded correctly. If not, the script is prevented from further running.
if [[ $(module list | grep 'intel/2023.2.1') == ""  || $(module list | grep 'python/3.11.5') == "" || $(module list | grep 'vasp/6.4.2') ]]; then
	module list
	echo "Your modules are not loaded correctly. Cancelling job... "
	exit 1
else
	module list
	echo "Your modules are loaded correctly. Proceeding to activate ASE..."
fi

# This loads the python virtualenv for VASP, using ASE
function load_ase() {
	source ~/software/python/virtualenvs/ase/bin/activate
}
load_ase

# Source .bashrc for environment variables and aliases
function bashrc() {
	source ~/.bashrc
}
bashrc

echo "################################
      # We are now Running your job! #
      ################################"
# Use of date variable allows the job to automatically find the right folder for the jobs at the time of submission
#python ~/scratch/submitted_date/$date/job1/job1.py
#python ~/scratch/submitted_date/$date/job2/job2.py
