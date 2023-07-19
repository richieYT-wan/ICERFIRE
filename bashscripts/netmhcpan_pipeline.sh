#! /bin/bash

# Assuming the input file is comma-separated, without header or index,
# with format
# peptide,wildtype,hla,target
# where target is optional

input_file=${1}
filename=$(basename ${input_file})
basenm="${filename%.*}"
final_fn="${basenm}_scored_output"

# TODO: HERE PATHS ARE TO BE REPLACED BY /tools/src/ICERFIRE-1.0/
OUTDIR=${2}
# TODO: HERE PATH IS TO BE REPLACED BY THE PATH TO NETMHCPAN4-1 (command line) in HealthTechCluster
NETMHCPAN=${3}
KERNDIST=${4}

echo "#######################"
echo "Processing ICOREs with NetMHCpan"
echo "#######################"

# Saving the Peptides to a newfile to be used later

echo "Peptide" > "${OUTDIR}base_file.txt"
# loop over lines in input file
line_number=1
while IFS=',' read -r column1 column2 column3; do
  ########################
  #   PROCESSING mutant  #
  ########################


  # create temporary file
  tmp_file_mut=$(mktemp)
  # write first column to temporary file
  echo ">seq$line_number" > "$tmp_file_mut"
  echo "$column1" >> "$tmp_file_mut"
  # Save the peps to the base file used to re-cover the full peptides
  echo "$column1" >> "${OUTDIR}base_file.txt"

  # run program with temporary file and second column as arguments
  ${NETMHCPAN} -l 8,9,10,11,12 -f "$tmp_file_mut" -a "$column3" -s > "${OUTDIR}output_tmp_mut.txt"

  # discard lines starting with #, select the line that has the lowest value for the Rnk_EL column, discard blank lines and lines starting with HLA, and write output to final_output.txt
  # Saving only selected columns
  # Pos HLA Peptide Core Of Gp Gl Ip Il Icore Rnk_EL
  grep -v '^#' "${OUTDIR}output_tmp_mut.txt" | grep -v '^HLA' | grep -v '^\-' | awk 'NF' | awk 'FNR==2{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$13}' >> "${OUTDIR}almostfinal_output_mut.txt"

  ########################
  # PROCESSING WILD TYPE #
  ########################

  # Get the starting position, and length to get the wt_aligned_icore
  start=$(grep -v '^#' "${OUTDIR}output_tmp_mut.txt" | grep -v '^HLA' | grep -v '^\-' | awk 'NF'| head -n 2 | awk '{print $1}' | awk 'FNR==2{print $1}')
  length=$(grep -v '^#' "${OUTDIR}output_tmp_mut.txt" | grep -v '^HLA' | grep -v '^\-' | awk 'NF'| head -n 2 | awk 'FNR==2{print length($10)}')
  wt_aligned="${column2:$((start-1)):length}"

  # Saving temp file and running netmhcpan in peptide mode
  tmp_file_wt=$(mktemp)
  echo "$wt_aligned" > "$tmp_file_wt"
  ${NETMHCPAN} -l ${length} -p "$tmp_file_wt" -a "$column3" -s > "${OUTDIR}output_tmp_wt.txt"

  # Discard line and grepping, extract only the aligned peptide and rank
  grep -v '^#' "${OUTDIR}output_tmp_wt.txt" | grep -v '^HLA' | grep -v '^\-' | awk 'NF'| awk 'FNR==2{print $3,$13}' >> "${OUTDIR}almostfinal_output_wt.txt"
  # remove temporary file
  rm "$tmp_file_mut"
  rm "$tmp_file_wt"

  # increment line number
  line_number=$((line_number+1))
done < "$input_file"


# Writing header and results to final output mut
echo "icore_start_pos HLA Pep Core Of Gp Gl Ip Il icore_mut EL_rank_mut"> "${OUTDIR}final_output_mut.txt"
cat "${OUTDIR}almostfinal_output_mut.txt" >> "${OUTDIR}final_output_mut.txt"

# Doing the same to final output WT
echo "icore_wt_aligned EL_rank_wt_aligned" > "${OUTDIR}final_output_wt.txt"
cat "${OUTDIR}almostfinal_output_wt.txt" >> "${OUTDIR}final_output_wt.txt"

# Pasting and saving final output
paste -d' ' "${OUTDIR}final_output_mut.txt" "${OUTDIR}final_output_wt.txt" > "${OUTDIR}merged_output.txt"

paste -d' ' "${OUTDIR}base_file.txt" "${OUTDIR}merged_output.txt" > "${OUTDIR}merged_final_output_tmp.txt" && mv "${OUTDIR}merged_final_output_tmp.txt" "${OUTDIR}merged_output.txt"

# Deleting all temp files
rm "${OUTDIR}almostfinal_output_wt.txt"
rm "${OUTDIR}almostfinal_output_mut.txt"
rm "${OUTDIR}final_output_wt.txt"
rm "${OUTDIR}final_output_mut.txt"
rm "${OUTDIR}output_tmp_mut.txt"
rm "${OUTDIR}output_tmp_wt.txt"
rm "${OUTDIR}base_file.txt"


echo " "
echo "#######################"
echo "Done processing ICOREs"
echo "#######################"

echo " "

echo "#######################"
echo "Processing Dissimilarity score"
echo "#######################"

echo "${KERNDIST}"
len1=$(wc -l "${OUTDIR}merged_output.txt" | awk '{print $1}')
fn=$(basename "${OUTDIR}merged_output.txt")
echo "icore_similarity_score" > "${OUTDIR}${fn}.kerndist"
for i in $(seq 2 $len1)
do
  awk -v line=$i 'FNR==line {print $11}' "${OUTDIR}merged_output.txt" > tmp1.pep
  awk -v line=$i 'FNR==line {print $13}' "${OUTDIR}merged_output.txt" > tmp2.pep
  ${KERNDIST} -blf "/tools/src/ICERFIRE-1.0/data/Matrices/blosum62.qij" -kmin 3 -kmax 8 ./tmp1.pep ./tmp2.pep | tail -n 3 >> "${OUTDIR}${fn}.kerndist"
  #{origin = $1} {$1 = $origin; print}
done
#rm tmp1.pep tmp2.pep
#paste -d' ' "${OUTDIR}merged_output.txt" "${OUTDIR}${fn}.kerndist" > "${OUTDIR}${final_fn}.txt"
#rm "${OUTDIR}merged_output.txt" "${OUTDIR}${fn}.kerndist"
#awk -F ' ' 'NR>1 {print $13}' "${OUTDIR}${final_fn}.txt" > "${OUTDIR}${final_fn}_wt_icore.txt"



