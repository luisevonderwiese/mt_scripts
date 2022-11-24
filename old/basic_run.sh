#!/bin/bash
FOLDER="../data/language_alignments"
#ALIGNMENT_FILE="Indo-European_WALS_BinaryOutgroup_March21_stability.ONLYFIN_common_grammatical.phy"
#ALIGNMENT_FILE="Indo-European_WALS_BinaryOutgroup_March21_stability.ONLYFIN_common_grammatical_filtered_-5.0.phy"



INPUT=${FOLDER}/${ALIGNMENT_FILE}
#DIR="../output/raxml/MORPHOLOGICAL/"
#DIR="../output/raxml/MORPHOLOGICAL_FILTERED/"
mkdir $DIR
PREFIX="${DIR}T_"
MODEL="BIN+G"
./../tools/raxml-ng/build/bin/raxml-ng --check --msa $INPUT --model $MODEL --prefix ${PREFIX}_1
./../tools/raxml-ng/build/bin/raxml-ng --parse --msa $INPUT --model $MODEL --prefix ${PREFIX}_2
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2
#./../raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX}_4 --threads 2 --seed 2 --tree pars{25},rand{25}
#./../tools/raxml-ng/build/bin/raxml-ng --search1 --msa $INPUT --model $MODEL --prefix ${PREFIX}_5 --threads 2 --seed 2
#grep "Final LogLikelihood:" ${PREFIX}_{3,4,5}.raxml.log
#./../tools/raxml-ng/build/bin/raxml-ng --rfdist --tree ${PREFIX}.raxml.mlTrees --prefix ${PREFIX}_RF
