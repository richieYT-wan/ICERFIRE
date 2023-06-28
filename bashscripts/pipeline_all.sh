#! /usr/bin/bash

source /home/projects/vaccine/people/yatwan/anaconda3/etc/profile.d/conda.sh
source activate pynn

FILENAME=${1}
filename=$(basename ${FILENAME})
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

echo " "
echo "#######################"
echo "Processing PepX score"
echo "#######################"
sh query_pepx.sh "${TMP}${final_fn}_wt_icore.txt"

cd ${PYDIR}
echo " "
echo "#######################"
echo " Running Model"
echo "#######################"
python3 run_model.py -f "${TMP}${final_fn}.txt" -pf "${TMP}${final_fn}_wt_icore_pepx_output.csv"
#rm ${TMP}/*