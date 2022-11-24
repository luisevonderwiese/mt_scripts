#!/bin/bash
cat  "../data/trees/cognate_ml.tree" "../data/trees/cognate_ie_compatible.tree"  > cognates.trees
./../raxml-ng/build/bin/raxml-ng --rfdist --tree cognates.trees --prefix temp --nofiles
rm cognates.trees

