#!/bin/bash

module purge
module load StdEnv/2023 intel/2023.2.1 openmpi/4.1.5
module load vasp/6.4.2
module load python/3.11.5 scipy-stack

if [[ $(module list | grep 'intel/2023.2.1') == ""  || $(module list | grep 'python/3.11.5') == "" || $(module list | grep 'vasp/6.4.2') == "" ]]; then
	echo "Your modules are not loaded correctly. Cancelling job... "
	exit 1
else
	echo "Your modules are loaded correctly. Proceeding to activate ASE..."
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
