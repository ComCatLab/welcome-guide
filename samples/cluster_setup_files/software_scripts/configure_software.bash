#!/usr/bin/env bash
# This script creates templated configuration files for ASE, ccu, and autojob.
# It is intended for use as an executable via the --software-script option of
# cluster-setup.
# Usage:
#   configure_software.bash SUPPORT CONFIG
#       where the SUPPORT and CONFIG arguments are the desired support file and
#       configuration home directories which default to the absolute paths
#       represented by ./support_file_home and ./.config, respectively


function template-file ()
{
    file=$1
    to_replace=$2
    replacement=$3

    sed -i "s|${to_replace}|${replacement}|g" "$file"
}

support_file_home=${1:-$(realpath .)/support_file_home}
config_home=${2:-$(realpath .)/.config}
support_file_placeholder='{{ support_file_home }}'
autojob_placeholder='{{ template_dir }}'
mkdir -p "$config_home"

# Configure ASE
ase_home="$config_home"/ase
mkdir "$ase_home"
# Change this relative path to point to the ASE configuration file template
ase_config_template=templates/configuration/ase.ini.j2
ase_config="${ase_home}/ase.ini"
cp -v "$ase_config_template" "$ase_config"
template-file "$ase_config" "$support_file_placeholder" "$support_file_home"
{ echo ""; echo "export ASE_CONFIG_PATH=$ase_config/ase.ini"; } >> ~/.bashrc

# Configure ccu
ccu_home="$config_home"/ccu
mkdir "$ccu_home"

# Configure autojob
autojob_home="$config_home"/autojob
autojob_template_dir="$autojob_home"/templates
autojob_templates=autojob_templates.tar.zst
# Change this relative path to point to the autojob template archive
autojob_template_archive=support_files/"$autojob_templates"
# Change this relative path to point to the autojob configuration file template
autojob_config_template=templates/configuration/autojob.toml.j2
autojob_config="$config_home"/autojob/autojob.toml
mkdir -p "$autojob_template_dir"
cp -v "$autojob_template_archive" "$autojob_template_dir"
cp -v "$autojob_config_template" "$autojob_config"
template-file "$autojob_config" "$autojob_placeholder" "${autojob_template_dir}"
cd "$autojob_template_dir" || exit
tar -vxf "$autojob_templates" && rm -v "$autojob_templates"
cd ~- || exit
