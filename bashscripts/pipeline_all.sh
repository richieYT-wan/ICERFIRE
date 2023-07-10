#! /bin/bash
# TODO: Remove these placeholders as I only used them to test my script

FILENAME=${1}
USER_EXPR=${2}
ADD_EXPR=${3}
filename=$(basename ${FILENAME})
basenm="${filename%.*}"
final_fn="${basenm}_scored_output"

# USERDIR=/home/projects/vaccine/people/yatwan/
USERDIR="/tools/src/"
BASHDIR="${USERDIR}ICERFIRE-1.0/bashscripts/"
#USERDIR=/Users/riwa/Documents/code/
DATADIR="${USERDIR}ICERFIRE-1.0/data/"
TMP="${USERDIR}ICERFIRE-1.0/tmp/"
#NETMHCPAN=/home/projects/vaccine/people/morni/netMHCpan-4.1/netMHCpan
NETMHCPAN="/tools/src/netMHCpan-4.1/netMHCpan"
KERNDIST="${USERDIR}ICERFIRE-1.0/bin/pep_kernel_dist"
PEPXDIR="/home/databases/userdb/pepx/"
PYTHON="/home/ctools/opt/anaconda3_202105/bin/python3"
PYDIR="${USERDIR}ICERFIRE-1.0/pyscripts/"

# Go to the bashdir and run the bash commands
cd ${BASHDIR}
bash netmhcpan_pipeline.sh ${FILENAME} ${TMP} ${NETMHCPAN} ${KERNDIST}

# TODO: USER_EXPR should come from a checkbox in the front end (asking whether the user is providing expr values, false by default)
# TODO: ADD_EXPR should come from a checkbox in the front end (Asking whether expression should be added to the model)
# TODO: ADD_EXPR should be checked / true by default in the front end.
if [ "$USER_EXPR" = "false" ] && [ "$ADD_EXPR" = "true" ]; then
  echo " "
  echo "#######################"
  echo "Processing PepX score"
  echo "#######################"
  bash query_pepx.sh "${TMP}${final_fn}_wt_icore.txt"
  PF="${TMP}${final_fn}_wt_icore_pepx_output.csv"
elif [ "$USER_EXPR" = "true" ]; then
  # TODO: Here merge user expr (4th column) to the file
  echo 'total_gene_tpm' > "${TMP}{final_fn}_tmp_expr.txt"
  awk -F ',' 'BEGIN{OFS=","} {print $4}' sample_data_expr.txt >> "${TMP}${final_fn}_expr.txt"
  paste -d ' ' "${TMP}${final_fn}.txt" "${TMP}${final_fn}_tmp_expr.txt" > "${TMP}${final_fn}_tmp_merged.txt" && mv "${TMP}${final_fn}_tmp_merged.txt" "${TMP}${final_fn}.txt"
elif [ "$ADD_EXPR" = "false" ];then
  echo "User-provided expression values or no expression added ; Skipping PepX query"
  PF="None"
fi

# Go to the Python dir and run the final model script
cd ${PYDIR}
echo "$PF"
echo "#######################"
echo " Running Model"
echo "#######################"
$PYTHON run_model.py -f "${TMP}${final_fn}.txt" -pf "$PF" -ae "$ADD_EXPR"
