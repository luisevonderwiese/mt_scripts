#!/bin/bash
FOLDER="data/language_alignments"
ALIGNMENT_NAME="morpho_filtered_pars_final"
ALIGNMENT_FILE=${ALIGNMENT_NAME}".phy"




INPUT=${FOLDER}/${ALIGNMENT_FILE}
mkdir ${FOLDER}/${ALIGNMENT_NAME}_pars
PREFIX=${FOLDER}/${ALIGNMENT_NAME}_pars/pars
MODEL="BIN+G"
./../../tools/raxml-ng/build/bin/raxml-ng  --start --tree pars{500} --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2
