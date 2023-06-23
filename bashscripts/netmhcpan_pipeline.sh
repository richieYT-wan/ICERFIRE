#! /usr/bin/bash
DATADIR=/home/projects/vaccine/people/yatwan/netmhcpan/score_pipeline/data/
OUTDIR=/home/projects/vaccine/people/yatwan/netmhcpan/score_pipeline/output/
NETMHCPAN=/home/projects/vaccine/people/morni/netMHCpan-4.1/netMHCpan
KERNDIST=/home/projects/vaccine/people/morni/bin/pep_kernel_dist
# Assuming the input file is comma-separated txt file, without header or index,
# with format
# peptide,wildtype,hla,target
# where target is optional

input_file="${DATADIR}$1"
final_fn=$(basename ${input_file})
echo "#######################"
echo "Processing ICOREs with NetMHCpan"
echo "#######################"

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

  # run program with temporary file and second column as arguments
  ${NETMHCPAN} -l 8,9,10,11,12 -f "$tmp_file_mut" -a "$column3" -s > output_tmp_mut.txt

  # discard lines starting with #, select the line that has the lowest value for the Rnk_EL column, discard blank lines and lines starting with HLA, and write output to final_output.txt
  # Saving only selected columns
  # Pos HLA Peptide Core Of Gp Gl Ip Il Icore Rnk_EL
  grep -v '^#' output_tmp_mut.txt | grep -v '^HLA' | grep -v '^\-' | awk 'NF' | awk 'FNR==2{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$13}' >> "${OUTDIR}almostfinal_output_mut.txt"

  ########################
  # PROCESSING WILD TYPE #
  ########################

  # Get the starting position, and length to get the wt_aligned_icore
  start=$(grep -v '^#' output_tmp_mut.txt | grep -v '^HLA' | grep -v '^\-' | awk 'NF'| head -n 2 | awk '{print $1}' | awk 'FNR==2{print $1}')
  length=$(grep -v '^#' output_tmp_mut.txt | grep -v '^HLA' | grep -v '^\-' | awk 'NF'| head -n 2 | awk 'FNR==2{print length($10)}')
  wt_aligned="${column2:$((start-1)):length}"

  # Saving temp file and running netmhcpan in peptide mode
  tmp_file_wt=$(mktemp)
  echo "$wt_aligned" > "$tmp_file_wt"
  ${NETMHCPAN} -l ${length} -p "$tmp_file_wt" -a "$column3" -s > output_tmp_wt.txt

  # Discard line and grepping, extract only the aligned peptide and rank
  grep -v '^#' output_tmp_wt.txt | grep -v '^HLA' | grep -v '^\-' | awk 'NF'| awk 'FNR==2{print $3,$13}' >> "${OUTDIR}almostfinal_output_wt.txt"
  # remove temporary file
  rm "$tmp_file_mut"
  rm "$tmp_file_wt"

  # increment line number
  line_number=$((line_number+1))
done < "$input_file"


# Writing header and results to final output mut
echo "icore_start_pos HLA Peptide Core Of Gp Gl Ip Il icore_mut EL_rank_mut"> "${OUTDIR}final_output_mut.txt"
cat "${OUTDIR}almostfinal_output_mut.txt" >> "${OUTDIR}final_output_mut.txt"

# Doing the same to final output WT
echo "icore_wt_aligned EL_rank_wt_aligned" > "${OUTDIR}final_output_wt.txt"
cat "${OUTDIR}almostfinal_output_wt.txt" >> "${OUTDIR}final_output_wt.txt"

# Pasting and saving final output
paste -d' ' "${OUTDIR}final_output_mut.txt" "${OUTDIR}final_output_wt.txt" > "${OUTDIR}merged_output.txt"

# Deleting all temp files
rm "${OUTDIR}almostfinal_output_wt.txt"
rm "${OUTDIR}almostfinal_output_mut.txt"
rm "${OUTDIR}final_output_wt.txt"
rm "${OUTDIR}final_output_mut.txt"
rm output_tmp_mut.txt
rm output_tmp_wt.txt


echo " "
echo "#######################"
echo "Done processing ICOREs"
echo "#######################"

echo " "

echo "#######################"
echo "Processing Dissimilarity score"
echo "#######################"


len1=$(wc -l "${OUTDIR}merged_output.txt" | awk '{print $1}')
fn=$(basename "${OUTDIR}merged_output.txt")
echo ${fn}
echo "icore_similarity_score" > "${OUTDIR}${fn}.kerndist"
for i in $(seq 2 $len1)
do
  awk -v line=$i 'FNR==line {print $10}' "${OUTDIR}merged_output.txt" > tmp1.pep
  awk -v line=$i 'FNR==line {print $12}' "${OUTDIR}merged_output.txt" > tmp2.pep
   ${KERNDIST} -kmin 3 -kmax 8 ./tmp1.pep ./tmp2.pep | tail -n 1 | awk '{print $3}' | awk 'NR == 1 {print $1}' >> "${OUTDIR}${fn}.kerndist"
  #{origin = $1} {$1 = $origin; print}
done
rm tmp1.pep tmp2.pep
paste -d' ' "${OUTDIR}merged_output.txt" "${OUTDIR}${fn}.kerndist" > "${OUTDIR}${final_fn}.txt"
rm "${OUTDIR}merged_output.txt" "${OUTDIR}${fn}.kerndist"


echo "#######################"
echo "  Done Preprocessing"
echo "  Running Model"
echo "#######################"


