import os
import pandas as pd
from ete3 import Tree


def read_num_states():
    lines = open("temp/max_states.csv", 'r').read().split("\n")[1:-1]
    num_states = {}
    for line in lines:
        data = line.split(",")
        num_states[data[0]] = int(data[1])
    return num_states

def eval_lh(tree, name, model):
    tree.write(outfile="temp.tree")
    if model.startswith("BIN"):
        alignment_name = "morph_alignments/bin/" + name
    else:
        alignment_name = "morph_alignments/multi/" + name
    os.system('./../../tools/raxml-ng/build/bin/raxml-ng --evaluate --msa ' + alignment_name +
            ' --threads 2 --model ' + model + ' --tree temp.tree --prefix foo --nofiles' +
              " --opt-branches off " + '> out.txt')
    l_file = open('out.txt', 'r')
    lines = l_file.readlines()
    lh = 1
    for line in lines:
        if(line.startswith('Final LogLikelihood:')):
            lh = float(line.split(" ")[2].strip())
    #if lh == 0:
    #    print(open("out.txt", "r").read())
    os.remove("out.txt")
    os.remove("temp.tree")
    return lh


def calculate_cross_lhs():
    morph_data_multistate = pd.read_parquet("training_data/morph_data_multistate.parquet")
    morph_data_binarized = pd.read_parquet("training_data/morph_data_binarized.parquet")
    num_states_dict = read_num_states()
    lh_file = open("temp/lhs.csv", "w+")
    lh_file.write("alignment,eval_tree,model,lh\n")
    for index, row in morph_data_multistate.iterrows():
        multitree = Tree(row["newick_eval"])
        multiname = row['verbose_name']
        print(multiname)
        num_states = num_states_dict[multiname]
        multimodel = "MULTI" + str(num_states) + "_GTR"
        binname = multiname.split('.')[0] + ".BIN.phy"
        bintree =  Tree(morph_data_binarized.loc[(morph_data_binarized['verbose_name'] == binname)].iloc[0]["newick_eval"])
        bamt_lh =  eval_lh(multitree, binname, "BIN")
        lh_file.write(binname + "," + multiname + ",BIN," + str(bamt_lh) + "\n")
        mabt_lh = eval_lh(bintree, multiname, multimodel)
        lh_file.write(multiname + "," + binname + "," + multimodel + "," + str(mabt_lh) + "\n")
calculate_cross_lhs()
