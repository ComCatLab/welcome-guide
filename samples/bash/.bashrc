# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

bind '"\e[A":history-search-backward'
bind '"\e[B":history-search-forward'

complete -d cd pushd rmdir

alias sl='ls'
alias ll='ls -lh --color=tty'
alias showq='squeue'
alias ag='ase gui'
alias open="xdg-open"
alias msub="new_submission"    # organized sbatch, new_submission is located in my bin folder
alias beluga="ssh -Y username@beluga.computecanada.ca"
alias cedar="ssh -Y username@cedar.computecanada.ca"
alias graham="ssh -Y username@graham.computecanada.ca"
alias jobres="seff"
alias activate_env='source $COMP_CHEM_ENV'

# Show CWD in shell prompt
export PS1="\[\033[36m\]\u\[\033[m\]@\[\033[32m\]\h:\[\033[33;1m\]\w\[\033[m\]\$ "
export LC_ALL=en_US.UTF-8

# Setting paths
export USERNAME="$HOME/projects/def-samiras/username"
export PATH="$USERNAME/bin:$PATH"
export COMP_CHEM_ENV="$USERNAME/software/comp-chem-env-py3.11/bin/activate"

# Configure modules
module purge
module use "$USERNAME/modules"
