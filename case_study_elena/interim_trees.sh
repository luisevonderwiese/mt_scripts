#!/bin/bash
INPUT="data/language_alignments/$1"
MODEL="BIN+G"
PREFIX="data/trees/$2pars100rand100bing"
./../../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree pars{100},rand{100}
MODEL="BIN"
PREFIX="data/trees/$2pars100rand100bin"
./../../tools/raxml-ng/build/bin/raxml-ng --msa $INPUT --model $MODEL --prefix ${PREFIX} --threads 2 --seed 2 --tree pars{100},rand{100}
