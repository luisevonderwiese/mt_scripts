#!/bin/bash
FOLDER="data/language_alignments"
ALIGNMENT_FILE="morpho.phy"
INPUT=${FOLDER}/${ALIGNMENT_FILE}
DIR="/home/luise/master_thesis/scripts/output/standard_raxml/"
mkdir $DIR

MODEL="BINGAMMA"
PREFIX="testbing"
./../tools/standard-RAxML-master/raxmlHPC-AVX  -p 12345 -m $MODEL -s $INPUT -n $PREFIX -w $DIR -# 2

MODEL="BINCAT"
PREFIX="testbin"
./../tools/standard-RAxML-master/raxmlHPC-AVX  -p 12345 -m $MODEL -s $INPUT -n $PREFIX -w $DIR -# 2
