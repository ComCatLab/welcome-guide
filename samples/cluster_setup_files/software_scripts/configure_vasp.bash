#!/usr/bin/env bash
# This script copies and extracts VASP support files, including the PAW
# pseudopotentials, a Python script used to execute VASP, and the VDW kernel.
# It is intended for use as an executable via the --software-script option of
# cluster-setup.
# Usage:
#   configure_software.bash SUPPORT
#       where the SUPPORT argument is the desired support file home directory
#       and defaults to support_file_home

support_file_home=${1:-support_file_home}
mkdir "$support_file_home"
vasp_archive=vasp.tar.zst
cp -v support_files/"$vasp_archive" "$support_file_home"
cd "$support_file_home" || exit
tar -vxf "$vasp_archive" && rm -v "$vasp_archive"
cd ~- || exit
