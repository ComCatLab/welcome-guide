#!/bin/bash
# This script simply initiates an interactive session on the cluster to avoid using the login node.
# The session has 4 cores, 12 GB memory, X11 forwarding, and a 1 hour time limit.
(cd ~/scratch/;salloc --x11 --time=01:00:00 --mem-per-cpu=3G --ntasks=4 --account=$SALLOC_ACCOUNT)
