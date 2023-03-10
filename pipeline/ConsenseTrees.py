import os
import pandas as pd
import shutil


def calculate_consense_trees(data_type, model):
    if model == "BIN":
        d = "parquets/" + data_type + "/" + model
    else:
        d = "parquets/" + data_type + "/MULTI/" + model
    with os.scandir(d) as it:
        for entry in it:
            if not entry.is_dir():
                continue
            dir_path = os.path.join(d, entry.name)
            df_path = os.path.join(dir_path, "raxmlng_tree_data.parquet")
            if not os.path.exists(df_path):
                print(entry.name)
                #shutil.rmtree(dir_path)
                continue
            tree_path = os.path.join(dir_path, "all_eval.trees")
            df = pd.read_parquet(df_path)
            tree_file = open(tree_path, "w+")
            for index, row in df.iterrows():
                tree_file.write(row["newick_eval"] + "\n")
            os.system("./raxml-ng --redo --consense --tree " + tree_path + " --prefix " + os.path.join(dir_path, "consense"))


calculate_consense_trees("morph", "MK")
