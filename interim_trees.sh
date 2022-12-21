#!/bin/bash
FOLDER="data/language_alignments"
ALIGNMENT_FILE="morpho_filtered_union_best_nj.phy"
INPUT=${FOLDER}/${ALIGNMENT_FILE}
DIR="output/raxml/interim/"
mkdir $DIR

MODEL="BIN+G"
#PREFIX="${DIR}test"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree rand{1}
PREFIX="${DIR}filtered_pars100rand100bing"
./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree pars{100},rand{100}
MODEL="BIN"
PREFIX="${DIR}filtered_pars100rand100bin"
./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree pars{100},rand{100}
#PREFIX="${DIR}start_cognate"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree ../data/trees/cognate_ie_compatible.tree
#PREFIX="${DIR}pars{100},rand{100}_bin"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree pars{100},rand{100}
#PREFIX="${DIR}cognate_ie_compatible_bin"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree ../data/trees/cognate_ie_compatible.tree
#PREFIX="${DIR}geo_science_scaled_interval=[0,0.3]_bin"
#./../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree ../data/trees/geo_science_scaled_interval=[0,0.3].tree
