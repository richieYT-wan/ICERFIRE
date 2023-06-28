#! /usr/bin/bash

FILENAME=${1}
input_file="${DATADIR}${1}"
filename=$(basename ${input_file})
basenm="${filename%.*}"
final_fn="${basenm}_scored_output"
USERDIR=/home/projects/vaccine/people/yatwan/
#USERDIR=/Users/riwa/Documents/code/
DATADIR="${USERDIR}ICERFIRE/data/"
TMP="${USERDIR}ICERFIRE/tmp/"
NETMHCPAN=/home/projects/vaccine/people/morni/netMHCpan-4.1/netMHCpan
KERNDIST=/home/projects/vaccine/people/morni/bin/pep_kernel_dist
PEPXDIR="${USERDIR}pepx/"
PYDIR="${USERDIR}ICERFIRE/pyscripts/"

sh netmhcpan_pipeline.sh ${FILENAME}
sh query_pepx.sh "${TMP}${final_fn}_wt_icore.txt"