#!/bin/bash

for i in {001..999}; do
    filename=dunnielex_multi/$i.phy
    prefix=output/$i
    ./raxml-ng --msa $filename --model BIN --prefix $prefix --threads 2 --seed 2
done

i=1000
filename=dunnielex_multi/$i.phy
prefix=output/$i
./raxml-ng --msa $filename --model BIN --prefix $prefix --threads 2 --seed 2

#./raxml-ng --msa dunnielex.catg --model BIN --prefix output/catg --threads 2 --seed 2 --prob-msa on

