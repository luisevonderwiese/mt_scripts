#!/bin/bash

for d in morph_parquets/multi/*/; do
 ./raxml-ng --consense --tree ${d}all_eval.trees --prefix ${d}consense
done





