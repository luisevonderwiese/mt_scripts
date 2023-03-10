import os
import pandas as pd
from ete3 import Tree


def read_num_states(data_type):
    lines = open("temp/" + data_type + "/max_states.csv", 'r').read().split("\n")[1:-1]
    num_states = {}
    for line in lines:
        data = line.split(",")
        num_states[data[0]] = int(data[1])
    return num_states

def eval_lh(data_type, tree, name, model):
    tree.write(outfile="temp.tree")
    if data_type == "morph":
        alignment_name = "alignments/morph/"
        if model.startswith("BIN"):
            alignment_name += "/bin/" + name.split('.')[0] + '.BIN.phy'
        else:
            alignment_name += "/multi/" + name
    if data_type == "lang":
        print("Alignment names for lang!")
    os.system('./../../tools/raxml-ng/build/bin/raxml-ng --evaluate --msa ' + alignment_name +
            ' --threads 2 --model ' + model + ' --tree temp.tree --prefix foo --nofiles' +
              " --opt-branches off " + '> out.txt')
    l_file = open('out.txt', 'r')
    lines = l_file.readlines()
    lh = 1
    for line in lines:
        if(line.startswith('Final LogLikelihood:')):
            lh = float(line.split(" ")[2].strip())
    if lh == 1:
        print(open("out.txt", "r").read())
    os.remove("out.txt")
    os.remove("temp.tree")
    return lh

# !! maybe needs to be adapted for language
def read_df_for_model(data_type, model):
    if data_type == "morph":
        d = "training_data/morph"
        if model == "BIN":
            df = pd.read_parquet(os.path.join(d, "binarized.parquet"))
            names = []
            for index, row in df.iterrows():
                names.append(row['verbose_name'].split('.')[0] + '.phy')
            df['verbose_name'] = names
        if model == "GTR":
            df = pd.read_parquet(os.path.join(d, "MULTI_GTR.parquet"))
        if model == "MK":
            df = pd.read_parquet(os.path.join(d, "full_MK.parquet"))
    elif data_type == "lang":
        print("Adapt training data paths for lang!")
    return df

def get_concrete_model(model, num_states):
    if model == "BIN":
        return model
    if model == "GTR":
        return "MULTI" + str(num_states) + "_GTR"
    if model == "MK":
        return "MULTI" + str(num_states) + "_MK"

def calculate_cross_lhs(data_type, model1, model2, outfile):
    df1 = read_df_for_model(data_type, model1)
    df2 = read_df_for_model(data_type, model2)
    num_states_dict = read_num_states(data_type)
    for i, row in df1.iterrows():
        name = row['verbose_name']
        if name not in num_states_dict:
            continue
        print(name)
        num_states = num_states_dict[name]
        concrete_model1 = get_concrete_model(model1, num_states)
        concrete_model2 = get_concrete_model(model1, num_states)
        df2_sub = df2.loc[(df2['verbose_name'] == name)]
        if (len(df2_sub) == 0):
            print(name + " not in second df!")
        tree1 = Tree(row["newick_eval"])
        tree2 = Tree(df2_sub.iloc[0]["newick_eval"])
        lh_t1_a2 = eval_lh(data_type, tree1, name, concrete_model2)
        outfile.write(name + "," + model1  + "," + model2 + "," + str(lh_t1_a2) + "\n")
        lh_t2_a1 = eval_lh(data_type, tree2, name, concrete_model1)
        outfile.write(name + "," + model2  + "," + model1 + "," + str(lh_t2_a1) + "\n")

def calculate_cross_all_cross_lhs(data_type):
    outfile = open("temp/" + data_type + "/lhs.csv", "w+")
    outfile.write("name,tree_model,eval_model,lh\n")
    calculate_cross_lhs(data_type, "BIN", "GTR", outfile)
    calculate_cross_lhs(data_type, "BIN", "MK", outfile)
    calculate_cross_lhs(data_type, "GTR", "MK", outfile)


calculate_cross_all_cross_lhs("morph")


