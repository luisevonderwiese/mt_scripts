import pandas as pd

from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import os
import json


' Code for Converting special cases, cannot be generalized'


def convert_json(path, data_set_id):
    with open(path, 'r') as f:
        data = json.load(f)
    language_dict = {}
    matrix = {}
    for language_id in data['Languages']:
        language_dict[int(language_id)] = data['Languages'][language_id]['Name'].replace(" ", "")
        matrix[int(language_id)] = ""
    data_sets = data['TypologicalDataSets']

    grids = data_sets[data_set_id]['TypologicalGrids']
    for grid_id in grids:
        features = grids[grid_id]['TypologicalFeatures']
        for feature_id in features:
            variants = features[feature_id]['TypologicalFeatureVariants']
            for variant_id in variants:
                d = {}
                if not 'RecordedValues' in variants[variant_id]:
                    continue
                values = variants[variant_id]['RecordedValues']
                for value_id in values:
                    entry = values[value_id]
                    d[entry['FkLanguageId']] = entry['Value']
                for language_id in language_dict:
                    if language_id in d:
                        if d[language_id]:
                            matrix[language_id] += '1'
                        else:
                            matrix[language_id] += '0'
                    else:
                        matrix[language_id] += '-'
    records = [SeqRecord(matrix[language_id], id=language_dict[language_id]) for language_id in matrix]
    align = MultipleSeqAlignment(records, annotations={}, column_annotations={})
    file_name = path.split("/")[-1].split(".")[0] + ".phy"
    with open(output_dir + file_name,"w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(align)
    #align = AlignIO.read(output_dir + file_name, "phylip-relaxed")


def convert_diacl():
    input_dir = "conversion_inputs/json/"
    output_dir = "alignment_uploads/"
    with os.scandir(input_dir) as it:
        for entry in it:
            if not entry.is_file():
                continue
            path = input_dir + entry.name
            data_set_id = entry.name.split('.')[0].split('-')[-1]
            convert_json(path, data_set_id)

#convert_diacl()


def convert_mouton_atlas():
    file_name = "datasets/raw/open/Mouton Atlas/Mouton_Atlas_Appendix2c.xlsx"
    dfs = pd.read_excel(file_name, dtype='object')
    columns = dfs.keys()
    d = {}
    taxa = []
    full_length = 0
    for column in columns:
        if column == "Lang_name.1":
            break
        if column.startswith("Unnamed"):
            continue
        if column == "Lang_name":
            for name in dfs[column]:
                d[name] = ""
                taxa.append(name)
            continue
        current_length = 0
        for i,v in enumerate(dfs[column]):
            if v != "-" and current_length == 0:
                current_length = len(v)
                full_length += current_length
            if v == "-":
                if current_length == 0:
                    j = 1
                    while(True):
                        w = dfs[column][j]
                        if w != "-" and current_length == 0:
                            current_length = len(w)
                            break
                        j = j + 1
                    full_length += current_length
                for k in range(current_length-1):
                    v = v + "-"
            d[taxa[i]] = d[taxa[i]] + v

    records = [SeqRecord(d[key], id=key.replace("(", "").replace(")", "").replace(" ", "_")) for key in d.keys()]
    align = MultipleSeqAlignment(records, annotations={}, column_annotations={})
    file_name = "mouton_atlas.ms.phy"
    with open(file_name,"w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(align)

#convert_mouton_atlas()




#------------------ DRAFT CODE -------------------------------------------

def df_from_diacl_file(path, data_set_id):
    with open(path, 'r') as f:
        data = json.load(f)
    language_dict = {}
    for language_id in data['Languages']:
        language_dict[int(language_id)] = data['Languages'][language_id]['Name'].replace(" ", "")

    language_id_col = []
    char_id_col = []
    value_col = []

    data_sets = data['TypologicalDataSets']
    grids = data_sets[data_set_id]['TypologicalGrids']
    for grid_id in grids:
        features = grids[grid_id]['TypologicalFeatures']
        for feature_id in features:
            variants = features[feature_id]['TypologicalFeatureVariants']
            for variant_id in variants:
                d = {}
                if not 'RecordedValues' in variants[variant_id]:
                    continue
                values = variants[variant_id]['RecordedValues']
                for value_id in values:
                    language_id = values[value_id]['FkLanguageId']
                    value = values[value_id]['Value']
                    if value:
                    language_id_col.append(language_dict[language_id])
                    char_id_col.append(variant_id)
                    value_col.append(value)

    df = pd.DataFrame()
    df["Language_ID"] = language_id_col
    df["Char_ID"] = char_id_col
    df["Value"] = value_col

    return df

def convert_diacl():
    d = "raw/DiACL"
    with os.scandir(d) as it:
        for entry in it:
            if not entry.is_file():
                continue
            path = os.path.join(d, entry.name)
            data_set_id = entry.name.split('.')[0].split('-')[-1]
            df = df_from_diacl_file(path, data_set_id)
            align = get_bin_align(df)
            print(align)
            file_name = os.path.join(output_dir, entry.name.split(".")[0] +  ".phy")
            with open(file_name,"w+") as f:
                writer = RelaxedPhylipWriter(f)
                writer.write_alignment(align)
