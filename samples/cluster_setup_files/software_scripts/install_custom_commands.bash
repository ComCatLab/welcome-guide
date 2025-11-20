#!/usr/bin/env bash

software_home=$1

if [[ $software_home = "" ]]; then echo "Error: No software home directory specified"; exit; fi
if ! test -e "$software_home"; then echo "Error: Software home directory $software_home does not exist"; exit; fi

# Change this relative path to point to the custom commands archive
cp -v sources/custom_commands.tar.zst "$software_home"
cd "$software_home" || exit
tar -vxf "custom_commands.tar.zst" && rm -v "custom_commands.tar.zst"
cd ~- || exit
{ echo ""; echo "module load custom-commands"; } >> ~/.bashrc
