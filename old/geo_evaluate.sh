#!/bin/bash
TREE="../output/geo.tree"
FOLDER="../data/language_alignments"
#ALIGNMENT_FILE="IE2011_RelaxedCovarion_AllSingletonsGeo.phy"
ALIGNMENT_FILE="Indo-European_WALS_BinaryOutgroup_March21_stability.ONLYFIN_common_grammatical.phy"

INPUT=${FOLDER}/${ALIGNMENT_FILE}
#DIR="../output/raxml/COGNATE/"
DIR="../output/raxml/temp/"
mkdir $DIR
PREFIX="${DIR}E"
MODEL="BIN"

./../raxml-ng/build/bin/raxml-ng --evaluate --msa  $INPUT --threads 2 --model $MODEL --tree $TREE --prefix ${PREFIX}_1

grep "Final LogLikelihood:" ${PREFIX}_1.raxml.log | cut -d ' ' -f3
rm -rf $DIR

