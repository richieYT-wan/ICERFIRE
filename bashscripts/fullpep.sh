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
  echo $line_number

  # create temporary file
  tmp_file_mut=$(mktemp)
  # write first column to temporary file
  echo "$column1" > "$tmp_file_mut"

  # run program with temporary file and second column as arguments
  ${NETMHCPAN} -l 8,9,10,11,12 -f "$tmp_file_mut" -a "$column3" -s -p > output_tmp_mut.txt

  # discard lines starting with #, select the line that has the lowest value for the Rnk_EL column, discard blank lines and lines starting with HLA, and write output to final_output.txt
  # Saving only selected columns
  # Pos HLA Peptide Core Of Gp Gl Ip Il Icore Rnk_EL
  grep -v '^#' output_tmp_mut.txt | grep -v '^HLA' | grep -v '^\-' | awk 'NF' | awk 'FNR==2{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$13}' >> "${OUTDIR}almostfinal_FULLPEP_output_mut.txt"

  # remove temporary file
  rm "$tmp_file_mut"

  # increment line number
  line_number=$((line_number+1))

done < "$input_file"


