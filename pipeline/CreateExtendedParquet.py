import os
import pandas as pd
import matplotlib.pyplot as plt
from plotly import graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from ete3 import Tree
import numpy as np



def rf_distance(t1, t2):
    rf, max_rf, common_leaves, parts_t1, parts_t2,discard_t1, discart_t2 = t1.robinson_foulds(t2, unrooted_trees = True)
    if max_rf == 0:
        print("?!")
        return 0
    return rf/max_rf

def read_consensus_trees(data_type, model):
    if model == "BIN":
        d = "parquets/" + data_type + "/BIN/"
    elif model == "GTR":
        d = "parquets/" + data_type + "/MULTI/GTR/"
    elif model == "MK":
        d = "parquets/" + data_type + "/MULTI/MK/"
    else:
        print(model + " does not exist!")
    consensus_trees = {}
    with os.scandir(d) as it:
        for entry in it:
            if not entry.is_dir():
                continue
            tree_path = os.path.join(d, os.path.join(entry.name, "consense.raxml.consensusTreeMR"))
            if not os.path.exists(tree_path):
                print("No consensus tree for " + model + " and " + entry.name)
                tree = Tree()
            else:
                tree = Tree(tree_path)
            #tree.resolve_polytomy(recursive=True)
            name = entry.name.split(".")[0] + ".phy"
            consensus_trees[name] = tree
    return consensus_trees

#Original model: Under this model a best tree was calculated with 100 tree searches
#Cross model: Now determine the likelihood of the best tree under this model
def cross_data_csv(data_type):
    out_file = open("temp/" + data_type + "/lhs/all.csv", "w+")
    out_file.write("name,original_model,cross_model,lh\n")
    lines = open("temp/" + data_type + "/lhs/BIN_GTR.csv", "r").read().split("\n")[1:-1]
    for line in lines:
        data = line.split(",")
        if data[0].endswith(".BIN.phy"):
            cross_model = "BIN"
            original_model="GTR"
            name = data[1]
        else:
            cross_model = "GTR"
            original_model="BIN"
            name = data[0]
        out_file.write(name + "," + original_model  + "," + cross_model + "," + data[3] + "\n")
    lines = open("temp/" + data_type + "/lhs/BIN_MK.csv", "r").read().split("\n")[1:-1]
    for line in lines:
        data = line.split(",")
        if data[0].endswith(".BIN.phy"):
            cross_model = "BIN"
            original_model="MK"
            name = data[1]
        else:
            cross_model = "MK"
            original_model="BIN"
            name = data[0]
        out_file.write(name + "," + original_model  + "," + cross_model + "," + data[3] + "\n")
    lines = open("temp/" + data_type + "/lhs/GTR_MK.csv", "r").read().split("\n")[1:-1]
    for line in lines:
        data = line.split(",")
        if data[2].endswith("GTR"):
            cross_model = "GTR"
            original_model="MK"
            name = data[0]
        else:
            cross_model = "MK"
            original_model="GTR"
            name = data[0]
        out_file.write(name + "," + original_model  + "," + cross_model + "," + data[3] + "\n")



#The column cross_lh_x contains the likelihood of the best tree found under the tree model evalut
def add_cross_data(df, data_type, eval_model, tree_model):
    lhs_df = pd.read_csv("temp/" + data_type + "/lhs/all.csv")
    lhs_df = lhs_df[lhs_df["original_model"] == tree_model]
    lhs_df = lhs_df[lhs_df["cross_model"] == eval_model]
    d = {}
    for idx, row in lhs_df.iterrows():
        d[row["name"]] = float(row["lh"])
    cross_lhs = []
    diffs = []
    for idx, row in df.iterrows():
        name = row['verbose_name']
        if name in d:
            cross_lh = d[name]
            eval_lh = row["llh_eval"]
            diff = cross_lh - eval_lh
            cross_lhs.append(cross_lh)
            diffs.append(diff)
        else:
            print("For " + name + " no cross evaluation with original model " + tree_model + " and cross model " +eval_model)
            cross_lhs.append(float("nan"))
            diffs.append(float("nan"))
    df["cross_llh_" + tree_model] = cross_lhs
    df["cross_diff_" + tree_model] = diffs
    return df


def merge_dfs(data_bin, data_gtr, data_mk):
    data_bin.columns = 'BIN_' + data_bin.columns.values
    data_gtr.columns = 'GTR_' + data_gtr.columns.values
    data_mk.columns = 'MK_' + data_mk.columns.values
    data_bin = data_bin.rename(columns={'BIN_verbose_name': 'verbose_name'})
    data_gtr = data_gtr.rename(columns={'GTR_verbose_name': 'verbose_name'})
    data_mk = data_mk.rename(columns={'MK_verbose_name': 'verbose_name'})
    df = pd.merge(data_bin, data_gtr, on='verbose_name', how='inner')
    df = pd.merge(df, data_mk, on='verbose_name', how='inner')
    df = df[(df.BIN_cross_llh_GTR != 1) & (df.BIN_cross_llh_MK != 1)
           & (df.GTR_cross_llh_BIN != 1) & (df.GTR_cross_llh_MK != 1)
           & (df.MK_cross_llh_BIN != 1) & (df.MK_cross_llh_GTR != 1)]
    return df


def add_rf_data(df, data_type):
    consensus_trees_bin = read_consensus_trees(data_type, "BIN")
    consensus_trees_gtr = read_consensus_trees(data_type, "GTR")
    consensus_trees_mk = read_consensus_trees(data_type, "MK")
    consensus_dist_bin_gtr = []
    consensus_dist_bin_mk = []
    consensus_dist_gtr_mk = []
    eval_dist_bin_gtr = []
    eval_dist_bin_mk = []
    eval_dist_gtr_mk = []


    for idx, row in df.iterrows():
        name = row["verbose_name"]
        c_tree_bin = consensus_trees_bin[name]
        c_tree_gtr = consensus_trees_gtr[name]
        c_tree_mk = consensus_trees_mk[name]
        e_tree_bin = Tree(row["BIN_newick_eval"])
        e_tree_gtr = Tree(row["GTR_newick_eval"])
        e_tree_mk = Tree(row["MK_newick_eval"])

        consensus_dist_bin_gtr.append(rf_distance(c_tree_bin, c_tree_gtr))
        consensus_dist_bin_mk.append(rf_distance(c_tree_bin, c_tree_mk))
        consensus_dist_gtr_mk.append(rf_distance(c_tree_gtr, c_tree_mk))
        eval_dist_bin_gtr.append(rf_distance(e_tree_bin, e_tree_gtr))
        eval_dist_bin_mk.append(rf_distance(e_tree_bin, e_tree_mk))
        eval_dist_gtr_mk.append(rf_distance(e_tree_gtr, e_tree_mk))

    df["consensus_dist_BIN_GTR"] = consensus_dist_bin_gtr
    df["consensus_dist_BIN_MK"] = consensus_dist_bin_mk
    df["consensus_dist_GTR_MK"] = consensus_dist_gtr_mk
    df["eval_dist_BIN_GTR"] = eval_dist_bin_gtr
    df["eval_dist_BIN_MK"] = eval_dist_bin_mk
    df["eval_dist_GTR_MK"] = eval_dist_gtr_mk

    return df

def add_aic_scores(df, data_type, model):
    lines = open("temp/" + data_type + "/aic.scores", "r").read().split("\n")[1:-1]
    aic_dict = {}
    caic_dict = {}
    bic_dict = {}
    for line in lines:
        data = line.split(',')
        if data_type == "morph":
            name = data[0].split(".")[0] + ".phy"
        if data_type == "lang":
            name = data[0].split(".")[0] + ".phy"
        cur_model = data[1]
        if cur_model != model:
            continue
        aic_dict[name]=float(data[2])
        caic_dict[name]=float(data[3])
        bic_dict[name]=float(data[4])
    aic_column = []
    caic_column = []
    bic_column = []
    for i, row in df.iterrows():
        name = row["verbose_name"]
        if name in aic_dict:
            aic_column.append(aic_dict[name])
        else:
            aic_column.append(0)
            print("No AIC score for model " + model + " and " + name)
        if name in caic_dict:
            caic_column.append(caic_dict[name])
        else:
            caic_column.append(0)
            print("No cAIC score for model " + model + " and " + name)
        if name in bic_dict:
            bic_column.append(bic_dict[name])
        else:
            bic_column.append(0)
            print("No BIC score for model " + model + " and " + name)
    df["AIC"] = aic_column
    df["cAIC"] = caic_column
    df["BIC"] = bic_column
    return df


def add_avg_col_states(df, data_type):
    avg_col_states = {}
    lines = open("temp/" + data_type + "/avg_col_states.csv", "r").read().split("\n")[1:-1]
    for line in lines:
        data = line.split(",")
        if data_type == "morph":
            name = data[0]
        if data_type == "lang":
            name = data[0]
        avg_col_states[name] = float(data[1])
    column = []
    for i, row in df.iterrows():
        column.append(avg_col_states[row["verbose_name"]])
    df["avg_col_states"] = column
    return df

def add_max_states(df, data_type):
    max_states = {}
    lines = open("temp/" + data_type + "/max_states.csv", "r").read().split("\n")[1:-1]
    for line in lines:
        data = line.split(",")
        if data_type == "morph":
            name = data[0]
        if data_type == "lang":
            name = data[0]
        max_states[name] = float(data[1])
    column = []
    for i, row in df.iterrows():
        column.append(max_states[row["verbose_name"]])
    df["max_states"] = column
    return df

def extend_df(data_type):
    data_gtr = pd.read_parquet("training_data/" + data_type + "/MULTI_GTR.parquet")
    data_mk = pd.read_parquet("training_data/" + data_type + "/full_MK.parquet")
    data_bin = pd.read_parquet("training_data/" + data_type + "/binarized.parquet")
    if data_type == "morph":
        #remove binary datasets
        data_mk = data_mk.groupby(data_mk.state_type).get_group("multistate")
        #unify names
        names = []
        for index, row in data_bin.iterrows():
            names.append(row['verbose_name'].split('.')[0] + '.phy')
        data_bin['verbose_name'] = names

    data_gtr = add_cross_data(data_gtr, data_type, "GTR", "BIN")
    data_gtr = add_cross_data(data_gtr, data_type, "GTR", "MK")
    data_mk = add_cross_data(data_mk, data_type, "MK", "BIN")
    data_mk = add_cross_data(data_mk, data_type, "MK", "GTR")
    data_bin = add_cross_data(data_bin, data_type, "BIN", "GTR")
    data_bin = add_cross_data(data_bin, data_type, "BIN", "MK")

    #data_gtr = add_aic_scores(data_gtr, data_type, "GTR")
    #data_mk = add_aic_scores(data_mk, data_type, "MK")
    #data_bin = add_aic_scores(data_bin, data_type, "BIN")

    df = merge_dfs(data_bin, data_gtr, data_mk)
    df = add_rf_data(df, data_type)
    df = add_avg_col_states(df, data_type)
    df = add_max_states(df, data_type)
    df.to_parquet("training_data/" + data_type + "/extended.parquet")

cross_data_csv("morph")
extend_df("morph")
#
