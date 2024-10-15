#! /usr/bin/bash

# This the main ICERFIRE-1.0 script. It acts as the full pipeline, doing the NetMHCpan, KernDist, PepX query, and Python script
# Yat-tsai Richie Wan, Jul 2023

###############################################################################
#               GENERAL SETTINGS: CUSTOMIZE TO YOUR SITE
###############################################################################

# Define the characters that can be used
characters="abcdefghijkmnopqrstuvwxyzABCDEFGHJKLMNOPQRSTUVWXYZ0123456789"
# Generate a random index between 0 and 61 (total number of characters)
index=$((RANDOM % 60))
# Get the character at the generated index
first_char="${characters:index:1}"
# Generate the remaining 4 characters as a combination of the defined characters
rest_chars=$(head /dev/urandom | tr -dc "$characters" | head -c 4)
# Combine the first and remaining characters

JOBID="${first_char}${rest_chars}"

while getopts ":f:a:u:" opt; do
  case ${opt} in
    f )
      FILENAME=$OPTARG
      ;;
    a )
      ADD_EXPR=$OPTARG
      ;;
    u )
      USER_EXPR=$OPTARG
      ;;
    \? )
      echo "Usage: $0 -f <file_path> -a <add_expr> (true/false) -u <user_exp> (true/false)"
      exit 1
      ;;
    : )
      echo "Invalid option: -$OPTARG requires an argument"
      exit 1
      ;;
  esac
done


# Shift the processed options so that $1, $2, etc. now refer to non-option arguments
shift $((OPTIND - 1))

# Access the values using the variable names
echo "File Path: $file_path"
echo "Add Expression: $add_expr"
echo "User Expression: $user_exp"

if [ -z "$FILENAME" ] || [ -z "$ADD_EXPR" ] || [ -z "$USER_EXPR" ]; then
  echo "Error: Missing required variables. Please provide values for file_path, add_expr, and user_exp."
  echo "Usage: $0 -f <file_path> -a <add_expr> (true/false) -u <user_exp> (true/false)"
  exit 1
fi

filename=$(basename ${FILENAME})
basenm="${filename%.*}"
final_fn="${basenm}_scored_output"

# other settings
PLATFORM="${UNIX}_${AR}"

# Replace USERDIR with the directory path containing ICERFIRE-1.0
USERDIR="/tools/src/"
# Replace PEPXDIR with the path to the folder containing the database file pepx-export.db
PEPXDIR="/home/databases/userdb/pepx/"
# Replace NETMHCPAN with the path to the binary of netMHCpan-4.1
NETMHCPAN="/tools/src/netMHCpan-4.1/netMHCpan"


BASHDIR="${USERDIR}ICERFIRE-1.0/bashscripts/"
SRCDIR="${USERDIR}ICERFIRE-1.0/src/"
DATADIR="${USERDIR}ICERFIRE-1.0/data/"
OUT="${USERDIR}ICERFIRE-1.0/output/${JOBID}/"

KERNDIST="${USERDIR}ICERFIRE-1.0/bin/pep_kernel_dist"


# Replace PYTHON a path to your installation of python3 
PYTHON="/home/ctools/opt/anaconda3_202105/bin/python3"
PYDIR="${USERDIR}ICERFIRE-1.0/pyscripts/"

mkdir -p ${OUT}
chmod 755 $OUT
# Go to the bashdir and run the bash commands
cd ${BASHDIR}
bash netmhcpan_pipeline.sh ${FILENAME} ${OUT} ${NETMHCPAN} ${KERNDIST}


case "$ADD_EXPR-$USER_EXPR" in
  "true-true")
#    echo "User-provided expression values; Skipping PepX query"
##    echo 'total_gene_tpm' > "${OUT}${final_fn}_tmp_expr.txt"
##    awk -F ',' '{print $4}' "${FILENAME}" >> "${OUT}${final_fn}_tmp_expr.txt"
##    paste -d ' ' "${OUT}${final_fn}.txt" "${OUT}${final_fn}_tmp_expr.txt" > "${OUT}${final_fn}_tmp_merged.txt" && mv "${OUT}${final_fn}_tmp_merged.txt" "${OUT}${final_fn}.txt"
##    ;;
#
#    if awk -F ',' 'NF>=4 && $4 ~ /^[0-9]*(\.[0-9]*)?$/ {exit 0} END {exit 1}' "${FILENAME}"; then
#        echo 'total_gene_tpm' > "${OUT}${final_fn}_tmp_expr.txt"
#        awk -F ',' '{print $4}' "${FILENAME}" >> "${OUT}${final_fn}_tmp_expr.txt"
#        paste -d ' ' "${OUT}${final_fn}.txt" "${OUT}${final_fn}_tmp_expr.txt" > "${OUT}${final_fn}_tmp_merged.txt" && mv "${OUT}${final_fn}_tmp_merged.txt" "${OUT}${final_fn}.txt"
#    else
#        echo "User-provided expression was selected but no fourth column found. Running expression database query"
#        echo " "
#        bash query_pepx.sh "${OUT}${final_fn}_wt_icore.txt" ${OUT}
#        PF="${OUT}${final_fn}_wt_icore_pepx_output.csv"
#    fi
#    ;;
    echo "User-provided expression values; Skipping PepX query"
    first_line=$(head -n 1 "${FILENAME}")

    # Use awk to count the number of fields (columns) and check if the fourth column is numeric
    num_columns=$(echo "$first_line" | awk -F ',' '{print NF}')
    fourth_column_is_numeric=$(echo "$first_line" | awk -F ',' '{if ($4 + 0 == $4) print "yes"; else print "no"}')

    # Perform actions based on whether there are 4 columns and the fourth column is numeric
    if [ "$num_columns" -eq 4 ] && [ "$fourth_column_is_numeric" = "yes" ]; then
        echo 'total_gene_tpm' > "${OUT}${final_fn}_tmp_expr.txt"
        awk -F ',' '{print $4}' "${FILENAME}" >> "${OUT}${final_fn}_tmp_expr.txt"
        paste -d ' ' "${OUT}${final_fn}.txt" "${OUT}${final_fn}_tmp_expr.txt" > "${OUT}${final_fn}_tmp_merged.txt" && mv "${OUT}${final_fn}_tmp_merged.txt" "${OUT}${final_fn}.txt"
    else
        echo "User-provided expression was selected but no valid fourth column found. Running expression database query"
        echo " "
        bash query_pepx.sh "${OUT}${final_fn}_wt_icore.txt" ${OUT}
        PF="${OUT}${final_fn}_wt_icore_pepx_output.csv"
    fi
  ;;
  "true-false")
    bash query_pepx.sh "${OUT}${final_fn}_wt_icore.txt" ${OUT}
    PF="${OUT}${final_fn}_wt_icore_pepx_output.csv"
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
cd ${SRCDIR}

chmod 755 "/home/locals/tools/src/ICERFIRE-1.0/src/"
$PYTHON run_model.py -j ${JOBID} -f "${OUT}${final_fn}.txt" -pf "$PF" -ae "$ADD_EXPR" -o "${OUT}" -ue "$USER_EXPR"
