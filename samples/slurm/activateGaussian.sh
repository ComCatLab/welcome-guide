#!/bin/bash

# Gaussian does not have any dependencies :angel:
module purge
module load gaussian/g16.c01
module load python/3.11.5 scipy-stack

if [[ $(module list | grep 'intel/2023.2.1') == ""  || $(module list | grep 'python/3.11.5') == "" || $(module list | grep 'vasp/6.4.2') == "" ]]; then
	echo "Your modules are not loaded correctly for Gaussian. Cancelling job... "
	exit 1
else
	echo "Your modules are loaded correctly for Gaussian. Proceeding to activate ASE..."
fi

echo "Changing directory to ~/software/python/virtualenvs/ase ..."
cd ~/software/python/virtualenvs/ase

function load_ase() {
	source ~/software/python/virtualenvs/ase/bin/activate
}

if [[ $(pwd | grep 'ase') == */software/python/virtualenvs/ase ]]; then
	pwd
	echo "You are in the right location! Activating ase..."
	load_ase
else
	echo "Please ensure you have the correct directory struture (~/software/python/virtualenvs/ase)..."
	echo "Exiting"
	exit 1
fi

#!/bin/bash

# Gaussian does not have any dependencies :angel:
module purge
module load gaussian/g16.c01
module load python/3.11.5 scipy-stack

if [[ $(module list | grep 'intel/2023.2.1') == ""  || $(module list | grep 'python/3.11.5') == "" || $(module list | grep 'vasp/6.4.2') == "" ]]; then
	echo "Your modules are not loaded correctly for Gaussian. Cancelling job... "
	exit 1
else
	echo "Your modules are loaded correctly for Gaussian. Proceeding to activate ASE..."
fi

echo "Changing directory to ~/software/python/virtualenvs/ase ..."
cd ~/software/python/virtualenvs/ase

function load_ase() {
	source ~/software/python/virtualenvs/ase/bin/activate
}

if [[ $(pwd | grep 'ase') == */software/python/virtualenvs/ase ]]; then
	pwd
	echo "You are in the right location! Activating ase..."
	load_ase
else
	echo "Please ensure you have the correct directory struture (~/software/python/virtualenvs/ase)..."
	echo "Exiting"
	exit 1
fi
