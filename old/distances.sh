#!/bin/bash
FOLDER="../data/language_alignments"
ALIGNMENT_FILE="IE2011_RelaxedCovarion_AllSingletonsGeo.phy"
#ALIGNMENT_FILE="Indo-European_WALS_BinaryOutgroup_March21_stability.ONLYFIN_common_grammatical.phy"

INPUT=${FOLDER}/${ALIGNMENT_FILE}
DIR="../output/raxml/"
#DIR="../output/raxml/MORPHOLOGICAL/"
mkdir $DIR
PREFIX="${DIR}"
MODEL="BIN"

#cat  "../output/scaled_geo.tree" "../output/cognate.tree" "../output/raxml/MORPHOLOGICAL/T_mltrees"  > ${PREFIX}scaled_alltrees
#./../raxml-ng/build/bin/raxml-ng --rfdist --tree ${PREFIX}scaled_alltrees --prefix ${PREFIX}scaled_alldistances

cat  "../output/scaled_geo.tree" "../output/cognate.tree" "../output/raxml/MORPHOLOGICAL_FILTERED/T_mltrees"  > ${PREFIX}scaled_filtered_alltrees
./../raxml-ng/build/bin/raxml-ng --rfdist --tree ${PREFIX}scaled_filtered_alltrees --prefix ${PREFIX}scaled_filtered_alldistances

cat  "../output/raxml/MORPHOLOGICAL_FILTERED/T_3.raxml.bestTree" "../output/raxml/MORPHOLOGICAL/T_3.raxml.bestTree"  > ${PREFIX}Trees_3
./../raxml-ng/build/bin/raxml-ng --rfdist --tree ${PREFIX}Trees_3 --prefix ${PREFIX}Distances_3

cat  "../output/raxml/MORPHOLOGICAL_FILTERED/T_4.raxml.bestTree" "../output/raxml/MORPHOLOGICAL/T_4.raxml.bestTree"  > ${PREFIX}Trees_4
./../raxml-ng/build/bin/raxml-ng --rfdist --tree ${PREFIX}Trees_4 --prefix ${PREFIX}Distances_4

cat  "../output/raxml/MORPHOLOGICAL_FILTERED/T_5.raxml.bestTree" "../output/raxml/MORPHOLOGICAL/T_5.raxml.bestTree"  > ${PREFIX}Trees_5
./../raxml-ng/build/bin/raxml-ng --rfdist --tree ${PREFIX}Trees_5 --prefix ${PREFIX}Distances_5

