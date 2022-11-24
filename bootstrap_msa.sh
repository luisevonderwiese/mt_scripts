#!/bin/bash
FOLDER="data/language_alignments"
ALIGNMENT_NAME="nj_subset_step10000000"
ALIGNMENT_FILE=${ALIGNMENT_NAME}".phy"




INPUT=${FOLDER}/${ALIGNMENT_FILE}
DIR="data/bootstrapping/"
mkdir $DIR
PREFIX=${FOLDER}/${ALIGNMENT_NAME}_bs/bs
MODEL="BIN+G"
./../tools/raxml-ng/build/bin/raxml-ng --bsmsa --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 -bs-trees 100
