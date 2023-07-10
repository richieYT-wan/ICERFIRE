#! /usr/bin/bash

FILENAME=${1}
filename=$(basename ${FILENAME})
basenm="${filename%.*}"
final_fn="${basenm}_scored_output"

# USERDIR=/home/projects/vaccine/people/yatwan/
USERDIR="/tools/src/ICERFIRE/"
BASHDIR=/home/projects/vaccine/people/yatwan/ICERFIRE/bashscripts/
#USERDIR=/Users/riwa/Documents/code/
DATADIR="${USERDIR}ICERFIRE/data/"
TMP="${USERDIR}ICERFIRE/tmp/PATHTONEWTMP"
#NETMHCPAN=/home/projects/vaccine/people/morni/netMHCpan-4.1/netMHCpan
NETMHCPAN="/tools/src/netMHCpan-4.1/netMHCpan"
KERNDIST="${USERDIR}pep_kernel_dist"
PEPXDIR="/home/databases/userdb/pepx/"
PYTHON="/home/ctools/opt/anaconda3_202105/bin/python"
PYDIR="${USERDIR}ICERFIRE/pyscripts/"

# TODO: Remove these placeholders as I only used them to test my script
USER_EXPR=${2}
ADD_EXPR=${3}
# Go to the bashdir and run the bash commands
cd ${BASHDIR}
sh netmhcpan_pipeline.sh ${FILENAME} ${TMP} ${NETMHCPAN} ${KERNDIST}

# TODO: USER_EXPR should come from a checkbox in the front end (asking whether the user is providing expr values, false by default)
# TODO: ADD_EXPR should come from a checkbox in the front end (Asking whether expression should be added to the model)
# TODO: ADD_EXPR should be checked / true by default in the front end.
if [ "$USER_EXPR" = false ] && [ "$ADD_EXPR" = true ]; then
  echo " "
  echo "#######################"
  echo "Processing PepX score"
  echo "#######################"
  sh query_pepx.sh "${TMP}${final_fn}_wt_icore.txt"
  PF="${TMP}${final_fn}_wt_icore_pepx_output.csv"
elif [ "$USER_EXPR" = true ]; then
  # TODO: Here merge user expr (4th column) to the file
  echo "xd"
elif [ "$ADD_EXPR" = false ];then
  echo "User-provided expression values or no expression added ; Skipping PepX query"
  PF="None"
fi

# Go to the Python dir and run the final model script
cd ${PYDIR}
echo " "
echo "#######################"
echo " Running Model"
echo "#######################"
PYTHON run_model.py -f "${TMP}${final_fn}.txt" -pf "$PF" -ae "$ADD_EXPR"
