import os
import pandas as pd
from ete3 import Tree


from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# All this code exists only because the taxa of BIN and MULTI alignments do not have the same names. This is shit, needs to be fixed and this code removed

def read_dict():
    lines = open("temp/lang/gerardi_and_reichert2021.BIN.cc.csv", "r").readlines()[1:]
    d={}
    for line in lines:
        data = line.split(",")
        if (len(data) != 2):
            print(line)
            continue
        d[data[1][:-1]] = data[0]
    return d

def rename_taxa_lang_bin(alignment_name):
    align = AlignIO.read(alignment_name, "phylip-relaxed")
    new_records = []
    for rec in align:
       new_records.append(SeqRecord(rec.seq, id=rec.id.split(".")[-1]))
    new_align = MultipleSeqAlignment(new_records, annotations={}, column_annotations={})
    new_alignemnt_path = "temp.phy"
    with open(new_alignemnt_path,"w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(new_align)

def gerardi_alignment(alignment_name):
    d = read_dict()
    align = AlignIO.read(alignment_name, "phylip-relaxed")
    new_records = []
    for rec in align:
       new_records.append(SeqRecord(rec.seq, id=d[rec.id]))
    new_align = MultipleSeqAlignment(new_records, annotations={}, column_annotations={})
    new_alignemnt_path = "temp.phy"
    with open(new_alignemnt_path,"w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(new_align)

def rename_taxa_lang_bin_in_tree(t):
    for leaf in t.iter_leaves():
        leaf.name = leaf.name.split(".")[-1]
    return t

def gerardi_tree(t):
    d = read_dict()
    for leaf in t.iter_leaves():
        leaf.name = d[leaf.name]
    return t

def get_neutral_lang_name(name):
    name_parts = name.split(".")
    del name_parts[1]
    return ".".join(name_parts)

def bin_name(neutral_name):
    name_parts = neutral_name.split(".")
    name_parts.insert(1, "BIN")
    return ".".join(name_parts)

def multi_name(neutral_name):
    name_parts = neutral_name.split(".")
    name_parts.insert(1, "MULTI")
    return ".".join(name_parts)

def add_neutral_names(df):
    neutral_names = []
    for i, row in df.iterrows():
        neutral_names.append(get_neutral_lang_name(row["verbose_name"]))
    df = df.rename(columns={'verbose_name': 'old_verbose_name'})
    df["verbose_name"] = neutral_names
    return df

def read_num_states(data_type):
    lines = open("temp/" + data_type + "/max_states.csv", 'r').read().split("\n")[1:-1]
    num_states = {}
    for line in lines:
        data = line.split(",")
        if data_type == "lang":
            num_states[get_neutral_lang_name(data[0])] = int(data[1])
        else:
            num_states[data[0]] = int(data[1])
    return num_states

def eval_lh(data_type, tree, name, model):
    tree.write(outfile="temp.tree")
    alignment_name = os.path.join("alignments", data_type)
    if data_type == "morph":
        if model.startswith("BIN"):
            alignment_name += "/bin/" + name.split('.')[0] + '.BIN.phy'
        else:
            alignment_name += "/multi/" + name
    if data_type == "lang":
        if model.startswith("BIN"):
            alignment_name += "/bin/" + bin_name(name)
            if name == "gerardi_and_reichert2021.cc.phy":
                gerardi_alignment(alignment_name)
            else:
                rename_taxa_lang_bin(alignment_name)
            alignment_name = "temp.phy"
        else:
            alignment_name += "/multi/" + multi_name(name)
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
    d = os.path.join("training_data", data_type)
    if data_type == "morph":
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
        df = pd.read_parquet(os.path.join(d, model + ".parquet"))
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
    if data_type == "lang":
        df1 = add_neutral_names(df1)
        df2 = add_neutral_names(df2)
    num_states_dict = read_num_states(data_type)
    for i, row in df1.iterrows():
        if model1 == "BIN" and model2 == "GTR" and i == 121:
            break
        name = row['verbose_name']
        if name not in num_states_dict:
            continue
        print(name)
        print(i)
        num_states = num_states_dict[name]
        concrete_model1 = get_concrete_model(model1, num_states)
        concrete_model2 = get_concrete_model(model2, num_states)
        df2_sub = df2.loc[(df2['verbose_name'] == name)]
        if (len(df2_sub) == 0):
            print(name + " not in second df!")
            continue
        tree1 = Tree(row["newick_eval"])
        if data_type == "lang" and model1 == "BIN":
            if name == "gerardi_and_reichert2021.cc.phy":
                tree1 =  gerardi_tree(tree1)
            else:
                tree1 =  rename_taxa_lang_bin_in_tree(tree1)
        tree2 = Tree(df2_sub.iloc[0]["newick_eval"])
        if data_type == "lang" and model2 == "BIN":
            if name == "gerardi_and_reichert2021.cc.phy":
                tree2 =  gerardi_tree(tree2)
            else:
                tree2 =  rename_taxa_lang_bin_in_tree(tree2)
        lh_t1_a2 = eval_lh(data_type, tree1, name, concrete_model2)
        outfile.write(name + "," + model1  + "," + model2 + "," + str(lh_t1_a2) + "\n")
        lh_t2_a1 = eval_lh(data_type, tree2, name, concrete_model1)
        outfile.write(name + "," + model2  + "," + model1 + "," + str(lh_t2_a1) + "\n")

def calculate_cross_all_cross_lhs(data_type):
    outfile = open("temp/" + data_type + "/lhs.csv", "w+")
    outfile.write("name,tree_model,eval_model,lh\n")
    calculate_cross_lhs(data_type, "BIN", "GTR", outfile)
    #calculate_cross_lhs(data_type, "BIN", "MK", outfile)
    #calculate_cross_lhs(data_type, "GTR", "MK", outfile)


calculate_cross_all_cross_lhs("lang")


