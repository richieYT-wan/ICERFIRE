#!/bin/bash

# This the main ICERFIRE-1.0 script. It acts as the full pipeline, doing the NetMHCpan, KernDist, PepX query, and Python script
# Yat, Jul 2023

###############################################################################
#               GENERAL SETTINGS: CUSTOMIZE TO YOUR SITE
###############################################################################

# determine where to store temporary files (must be writable to all users)

# Or maybe override
#export TMPDIR=/scratch

# determine platform
UNIX="Linux"
AR="x86_64"

# determine platform, do it through 'uname'
#UNIX=`uname -s`
#AR=`uname -m`

# WWWROOT of web server
WWWROOT=/var/www/html

# WWWpath to service
SERVICEPATH=/services/ICERFIRE-1.0

# other settings
PLATFORM="${UNIX}_${AR}"

if [ -z "$TMP" ]; then
	export TMP=/scratch
fi
#[ "$USER_EXPR" = "false" ] && [ "$ADD_EXPR" = "true" ]; th
# Expanding on blastdb
options=()
while (( $# > 0 )); do
   case $1 in
     "--jobid")
       shift
       JOBID=$1
     ;;
     "--add_expr")
      shift
      if [[ $1 == "yes" ]]; then
         ADD_EXPR="true"
      else
         ADD_EXPR="false"
      fi
      ;;
     "--usr_expr")
      shift
      if [[ $1 == "yes" ]]; then
         USER_EXPR="true"
      else
         USER_EXPR="false"
      fi
      ;;
     "--infile")
      shift
      FILENAME=$1
      ;;
#      *)
#       options+=("$1")
#      ;;
   esac
   shift
done

#options+=("-o")
#options+=("${WWWROOT}${SERVICEPATH}/tmp/${JOBID}")


FILENAME=${1}
#USER_EXPR=${2}
#ADD_EXPR=${3}
filename=$(basename ${FILENAME})
basenm="${filename%.*}"
final_fn="${basenm}_scored_output"


mkdir -p ${WWWROOT}${SERVICEPATH}/tmp/${JOBID}
mkdir -p /tmp/${JOBID}

USERDIR="/tools/src/"
BASHDIR="${USERDIR}ICERFIRE-1.0/bashscripts/"
DATADIR="${USERDIR}ICERFIRE-1.0/data/"
#TMP="${USERDIR}ICERFIRE-1.0/tmp/"
NETMHCPAN="/tools/src/netMHCpan-4.1/netMHCpan"
KERNDIST="${USERDIR}ICERFIRE-1.0/bin/pep_kernel_dist"
PEPXDIR="/home/databases/userdb/pepx/"
PYTHON="/home/ctools/opt/anaconda3_202105/bin/python3"
PYDIR="${USERDIR}ICERFIRE-1.0/pyscripts/"

# Go to the bashdir and run the bash commands
cd ${BASHDIR}
bash netmhcpan_pipeline.sh ${FILENAME} ${TMPDIR} ${NETMHCPAN} ${KERNDIST}

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
  echo "User-provided expression values; Skipping PepX query"
  echo 'total_gene_tpm' > "${TMP}${final_fn}_tmp_expr.txt"
  awk -F ',' '{print $4}' ${FILENAME} >> "${TMP}${final_fn}_tmp_expr.txt"
  paste -d ' ' "${TMP}${final_fn}.txt" "${TMP}${final_fn}_tmp_expr.txt" > "${TMP}${final_fn}_tmp_merged.txt" && mv "${TMP}${final_fn}_tmp_merged.txt" "${TMP}${final_fn}.txt"
elif [ "$ADD_EXPR" = "false" ];then
  echo "No expression added ; Skipping PepX query"
  PF="None"
fi

# Go to the Python dir and run the final model script
cd ${PYDIR}
echo "$PF"
echo "#######################"
echo " Running Model"
echo "#######################"
$PYTHON run_model.py -f "${TMP}${final_fn}.txt" -pf "$PF" -ae "$ADD_EXPR" -o "${WWWROOT}${SERVICEPATH}/tmp/${JOBID}"
rm "${TMP}*scored_output*"