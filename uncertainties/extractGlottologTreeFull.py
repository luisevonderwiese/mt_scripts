from ete3 import Tree
import re
import pandas as pd
import numpy as np


#%%

with open('../data/tree_glottolog_newick.txt') as f:
    raw = f.readlines()

trees = []

for i, ln in enumerate(raw):
    ln = ln.strip()
    ln = re.sub(r"\'[A-Z][^[]*\[", "[", ln)
    ln = re.sub(r"\][^']*\'", "]", ln)
    ln = re.sub(r"\[|\]", "", ln)
    ln = ln.replace(":1", "")
    trees.append(Tree(ln, format=1))

#%%


glot = Tree()
for t in trees:
    glot.add_child(t)


nonLeaves = [nd.name for nd in glot.traverse()
             if nd.name != '' and not nd.is_leaf()]

for i, nm in enumerate(nonLeaves):
    if i % 100 == 0:
        print(i)
    nd = glot & nm
    nd.name = ''
    nd.add_child(name=nm)

#%%

fn = "../data/glottolog-glottolog-cldf-ac0d616/cldf/languages.csv"
languages = pd.read_csv(fn)

#%%

gTaxa = languages.Glottocode[languages.ISO639P3code.notna()].values

tTaxa = glot.get_leaf_names()

gTaxa = np.intersect1d(gTaxa, tTaxa)

glot.prune([glot&x for x in gTaxa])


#%%

glot.write(outfile="glottolog.tre", format=9)

#%%
