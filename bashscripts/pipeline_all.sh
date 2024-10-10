#!/bin/bash

# This the main ICERFIRE-1.0 script. It acts as the full pipeline, doing the NetMHCpan, KernDist, PepX query, and Python script
# Yat, Jul 2023

###############################################################################
#               GENERAL SETTINGS: CUSTOMIZE TO YOUR SITE
###############################################################################

if [ -z "$TMP" ]; then
	export TMP=/scratch
fi


ADD_EXPR="false"
USER_EXPR="false"
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
       fi
       ;;
     "--usr_expr")
       shift
       if [[ $1 == "yes" ]]; then
         USER_EXPR="true"
       fi
       ;;
     "--infile")
       shift
       FILENAME=$1
       ;;
   esac
   shift
done


## TODO THIS UPDATE IS FOR COMMAND-LINE TESTING PURPOSES ONLY, UPDATE ACCORDINGLY WHEN DEPLOYING FOR THE SERVER
#FILENAME=${1}
#ADD_EXPR=${2}
#USER_EXPR=${3}

filename=$(basename ${FILENAME})
basenm="${filename%.*}"
final_fn="${basenm}_scored_output"

# determine platform
UNIX="Linux"
AR="x86_64"

# WWWROOT of web server
WWWROOT=/var/www/html

# WWWpath to service
SERVICEPATH=/services/ICERFIRE-1.0

# other settings
PLATFORM="${UNIX}_${AR}"
USERDIR="/tools/src/"
BASHDIR="${USERDIR}ICERFIRE-1.0/bashscripts/"
SRCDIR="${USERDIR}ICERFIRE-1.0/src/"
DATADIR="${USERDIR}ICERFIRE-1.0/data/"
TMP=${WWWROOT}${SERVICEPATH}/tmp/${JOBID}/
chmod 755 $TMP
#TMP="${USERDIR}ICERFIRE-1.0/tmp/"
NETMHCPAN="/tools/src/netMHCpan-4.1/netMHCpan"
KERNDIST="${USERDIR}ICERFIRE-1.0/bin/pep_kernel_dist"
PEPXDIR="/home/databases/userdb/pepx/"
PYTHON="/home/ctools/opt/anaconda3_202105/bin/python3"
PYDIR="${USERDIR}ICERFIRE-1.0/pyscripts/"

mkdir -p ${TMP}
mkdir -p /tmp/${JOBID}
# Go to the bashdir and run the bash commands

cat ${FILENAME} > ${TMP}RAW_INPUT_FILE.txt
cd ${BASHDIR}
bash netmhcpan_pipeline.sh ${FILENAME} ${TMP} ${NETMHCPAN} ${KERNDIST}


case "$ADD_EXPR-$USER_EXPR" in
  "true-true")
    echo "User-provided expression values; Skipping PepX query"
    first_line=$(head -n 1 "${FILENAME}")

    # Use awk to count the number of fields (columns) and check if the fourth column is numeric
    num_columns=$(echo "$first_line" | awk -F ',' '{print NF}')
    fourth_column_is_numeric=$(echo "$first_line" | awk -F ',' '{if ($4 + 0 == $4) print "yes"; else print "no"}')

    # Perform actions based on whether there are 4 columns and the fourth column is numeric
    if [ "$num_columns" -eq 4 ] && [ "$fourth_column_is_numeric" = "yes" ]; then
        echo 'total_gene_tpm' > "${TMP}${final_fn}_tmp_expr.txt"
        awk -F ',' '{print $4}' "${FILENAME}" >> "${TMP}${final_fn}_tmp_expr.txt"
        paste -d ' ' "${TMP}${final_fn}.txt" "${TMP}${final_fn}_tmp_expr.txt" > "${TMP}${final_fn}_tmp_merged.txt" && mv "${TMP}${final_fn}_tmp_merged.txt" "${TMP}${final_fn}.txt"
    else
        echo "User-provided expression was selected but no valid fourth column found. Running expression database query"
        echo " "
        bash query_pepx.sh "${TMP}${final_fn}_wt_icore.txt" ${TMP}
        PF="${TMP}${final_fn}_wt_icore_pepx_output.csv"
    fi
    ;;

  "true-false")
#    echo " "
#    echo "#######################"
#    echo "Processing PepX score"
#    echo "#######################"
    bash query_pepx.sh "${TMP}${final_fn}_wt_icore.txt" ${TMP}
    PF="${TMP}${final_fn}_wt_icore_pepx_output.csv"
    ;;

  "false-false")
    echo "No expression added; Skipping PepX query"
    PF="None"
    ;;

  "false-true")
    echo "Warning: Invalid combination - ADD_EXPR is 'false' but USER_EXPR is 'true'. Ignoring this part and continuing with the rest of the code."
    ;;

  *)
    echo "Invalid combination of flags: ADD_EXPR=$ADD_EXPR, USER_EXPR=$USER_EXPR"
    exit 1
    ;;
esac

# Go to the Python dir and run the final model script
#cd ${PYDIR}
cd ${SRCDIR}
#echo "HERE IS THE PF FILE $PF"
#echo " "
#echo "#######################"
#echo "     Running Model"
#echo "#######################"
chmod 755 "/home/locals/tools/src/ICERFIRE-1.0/src/"
$PYTHON run_model.py -j ${JOBID} -f "${TMP}${final_fn}.txt" -pf "$PF" -ae "$ADD_EXPR" -o "${TMP}" -ue "$USER_EXPR"
# > "${TMP}pylogs" 2>&1