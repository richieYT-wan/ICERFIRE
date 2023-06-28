#! /usr/bin/bash
#set -x
# sqlite> select * from expression_dataset where source like '%TCGA%' and title like '%PANCAN%';
# 1|TCGA|cancer-type|PANCAN|TCGA-PANCAN RNA expression|1|0| 
# TCGA pancan's dataset_id is 1 so here dataset_id=1 by default ; Can / will change behaviour if we start using other datasets

# Check if a filename argument is provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 <filepath> <dataset_id>"
  exit 1
fi

PEPXDIR='/home/projects/vaccine/people/yatwan/pepx/'
TMPDIR='/home/projects/vaccine/people/yatwan/ICERFIRE/tmp/'
filepath=$1
dataset_id=1
filename=$(basename "$filepath")
basenm="${filename%.*}"
database="${PEPXDIR}pepx-export.db"
output_file="${basenm}_pepx_output"

peptide_string=`cat $filepath | perl -ne 'while(<>){chomp; push @all, $_;} print join("\",\"",@all)'`

# Select only those fields, from the peptide_gene_tpms_collapsed table where dataset_id = 1 (id_1 == TCGA_PANCAN)
query='select dataset_id, peptide, total_peptide_tpm, total_scaled_peptide_tpm, total_gene_tpm, gene_symbols, gene_ensg_ids,'
query="$query gene_num_proteins, gene_num_matching_proteins, gene_frac_matching_proteins"
query="$query from peptide_gene_tpms_collapsed where"
query="$query dataset_id = $dataset_id"
query="$query and peptide in (\"$peptide_string\")"
query="$query order by peptide asc;"

echo "Running PepX query on ${1}, from ${database}, saving at ${output_file}"
sqlite3 $database -header "$query" > "${TMPDIR}${output_file}.csv"
echo "Query done ; Updating table format and moving temporary files"
# Replace | with commas to make it csv, using a temp file then mv to overwrite
sed 's/|/,/g' < "${TMPDIR}${output_file}.csv" > "${TMPDIR}${output_file}_temp.csv"
mv "${TMPDIR}${output_file}_temp.csv" "${TMPDIR}${output_file}.csv"
