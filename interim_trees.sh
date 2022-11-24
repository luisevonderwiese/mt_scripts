#!/bin/bash
FOLDER="data/language_alignments"
ALIGNMENT_FILE="morpho.phy"
INPUT=${FOLDER}/${ALIGNMENT_FILE}
DIR="../output/raxml/MORPHOLOGICAL_INTERIM/"
mkdir $DIR

MODEL="BIN+G"
PREFIX="${DIR}default_run"
./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2
#PREFIX="${DIR}pars_100_rand_100"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree pars{100},rand{100}
#PREFIX="${DIR}start_cognate"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree ../data/trees/cognate_ie_compatible.tree
#PREFIX="${DIR}pars{100},rand{100}_bin"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree pars{100},rand{100}
#PREFIX="${DIR}cognate_ie_compatible_bin"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree ../data/trees/cognate_ie_compatible.tree
#PREFIX="${DIR}geo_science_scaled_interval=[0,0.3]_bin"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree ../data/trees/geo_science_scaled_interval=[0,0.3].tree
