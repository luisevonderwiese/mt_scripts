import os
import pandas as pd
import math
import numpy as np
from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt
import json


states = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "!", "\"", "#", "$", "%", "&", "'", "(", ")", "*", "+",
",", "/", ":", ";", "<", "=", ">", "@", "[", "\\", "]", "^", "_", "{", "|", "}", "~"]
print(len(states))


input_dir = "raw/cldf"
output_dir = "alignment_uploads/alignments/"



def drop_unnecessary_columns(df):
    relevant_columns = ["ID", "Form_ID", "Language_ID", "Cognateset_ID", "Name", "Parameter_ID"]
    new_df = df
    for column in df.columns:
        if not column in relevant_columns:
            new_df = new_df.drop(column, axis=1)
    return new_df


def get_bin_align(df):
    languages = list(set(df['Language_ID']))
    print(str(len(languages)) + " languages")
    matrix = {}
    for language in languages:
        matrix[language] = ""
    chars = list(set(df['Char_ID']))
    chars.sort()
    num_chars = len(chars)
    print(str(num_chars) + " chars")
    for (c, char) in enumerate(chars):
        sub_df = df[df["Char_ID"] == char]
        char_languages = list(set(sub_df['Language_ID']))
        values = list(set(sub_df['Value']))
        values.sort()
        for language in languages:
            if language not in char_languages:
                bs = ['-' for _ in range(len(values))]
                matrix[language] += ''.join(bs)
                continue
            bs = ['0' for _ in range(len(values))]
            for(i, value) in enumerate(values):
                sub_sub_df = sub_df[(sub_df["Language_ID"] == language) & (sub_df["Value"] == value)]
                if(len(sub_sub_df.index) > 0):
                    bs[i] = '1'
            matrix[language] += ''.join(bs)
        print("(" + str(c+1) + "/" + str(num_chars)+")")
    records = [SeqRecord(matrix[language_id], id=str(language_id)) for language_id in matrix.keys()]
    align = MultipleSeqAlignment(records, annotations={}, column_annotations={})
    return align



def get_multi_align(df, drop_threshold, column_threshold):
    languages = list(set(df['Language_ID']))
    print(str(len(languages)) + " languages")
    matrix = {}
    for language in languages:
        matrix[language] = ""
    chars = list(set(df['Char_ID']))
    chars.sort()
    r = {
        "num_languages" : len(languages),
        "num_chars" : len(chars),
        "multi_cells" : 0,
        "multi_chars" : 0,
        "high_multi_chars" : 0,
        "high_state_chars" : 0,
        "dropped_chars" : 0,
        "max_states" : 0,
        "converted" : 0}
    for (c, char) in enumerate(chars):
        char_values = {}
        possible = True
        sub_df = df[df["Char_ID"] == char]
        char_languages = list(set(sub_df['Language_ID']))
        classes = list(set(sub_df['Value']))
        if (len(classes) > len(states)):
            r["high_state_chars"] += 1
            r["dropped_chars"] += 1
            continue
        counts = sub_df['Value'].value_counts()
        classes = sorted(classes, key=counts.get, reverse = True)
        column_multi_cells = 0
        for language in languages:
            sub_sub_df = sub_df[(sub_df["Language_ID"] == language)]
            if language not in char_languages or len(sub_sub_df) == 0:
                char_values[language] = '-'
                continue
            max_count = 0
            if len(sub_sub_df) > 1:
                column_multi_cells += 1
            for i, row in sub_sub_df.iterrows():
                if(counts[row['Value']] > max_count):
                    char_values[language] =  states[classes.index(row["Value"])]
            assert(language in char_values)
        if column_multi_cells / len(languages) > column_threshold:
            r["high_multi_chars"] += 1
            r["dropped_chars"] += 1
        else:
            for language in languages:
                matrix[language] += char_values[language]
        r["max_states"] = max(r["max_states"], len(classes))
        r["multi_cells"] += column_multi_cells
        if column_multi_cells > 0:
            r["multi_chars"] += 1
    if (r["dropped_chars"]/r["num_chars"] > drop_threshold):
        print(r)
        return (r, None)
    r["converted"] = 1
    print(r)
    records = [SeqRecord(matrix[language_id], id=str(language_id)) for language_id in matrix.keys()]
    align = MultipleSeqAlignment(records, annotations={}, column_annotations={})
    return (r, align)


def write_bin_catg(df, outfile):
    languages = list(set(df['Language_ID']))
    print(str(len(languages)) + " languages")
    matrix = {}
    chars = list(set(df['Char_ID']))
    chars.sort()
    tempfile = open("temp.txt", "w+")
    tempfile.write(" ".join(languages))
    tempfile.write("\n")
    num_sites = 0
    for (c, char) in enumerate(chars):
        value_counts = {}
        for language in languages:
            value_counts[language] = 0
        char_df = df[df["Char_ID"] == char]
        values = list(set(char_df['Value']))
        values.sort()
        #determine number of values of this characteristic per language
        for(i, value) in enumerate(values):
            value_df = char_df[char_df["Value"] == value]
            value_languages = list(set(value_df['Language_ID']))
            for language in languages:
                if language in value_languages:
                    value_counts[language] += 1
        #determine most probable cols and probablities and write to file
        num_sites += len(values)
        for(i, value) in enumerate(values):
            value_df = char_df[char_df["Value"] == value]
            value_languages = list(set(value_df['Language_ID']))
            col = ""
            probs = []
            for language in languages:
                if language in value_languages:
                    # only in this case 1 is more probable than 0
                    if value_counts[language] == 1:
                        col += "1"
                    else:
                        col += "0"
                    one_prob = 1 / value_counts[language]
                else:
                    col += "0"
                    # missing data
                    if value_counts[language] == 0:
                        probs.append("1.0,1.0")
                        continue
                        #alternative:
                        #one_prob = 1 / len(values)
                    else:
                        one_prob = 0.0
                zero_prob = 1.0 - one_prob
                one_prob = round(one_prob, 3)
                zero_prob = round(zero_prob, 3)
                probs.append(str(zero_prob)+","+str(one_prob))
            tempfile.write(col + " ")
            tempfile.write(" ".join(probs))
            tempfile.write("\n")
    outfile.write(str(len(languages)) + " ")
    outfile.write(str(num_sites) + "\n")
    tempfile.close()
    outfile.write(open("temp.txt", "r").read())
    os.remove("temp.txt")


def get_cognate_df(name):
    forms_df = pd.read_csv(os.path.join(input_dir, name + "/cldf/forms.csv"), low_memory=False)
    cognates_df = pd.read_csv(os.path.join(input_dir, name + "/cldf/cognates.csv"), low_memory=False)
    cognates_df = cognates_df[cognates_df.Cognate_Detection_Method == "expert"]
    forms_df = drop_unnecessary_columns(forms_df)
    cognates_df = drop_unnecessary_columns(cognates_df)

    full_df = pd.merge(forms_df, cognates_df, how="outer", left_on=['ID'], right_on=['Form_ID'])
    full_df = full_df[full_df.Cognateset_ID == full_df.Cognateset_ID]
    full_df = full_df.rename(columns={'Cognateset_ID': 'Value', 'Parameter_ID': 'Char_ID'})
    return full_df

def get_df_gerhard_csv(path):
    df = pd.read_csv(path)
    df = drop_unnecessary_columns(df)
    df = df.rename(columns={'Cognateset_ID': 'Value', 'Parameter_ID': 'Char_ID'})
    return df


def convert_cognate(name):
    df = get_cognate_df(name)
    align = get_bin_align(df)
    file_name = os.path.join(output_dir, name+ ".phy")
    with open(file_name,"w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(align)

def convert_cognate_multivalue(name, drop_threshold, column_threshold):
    df = get_cognate_df(name)
    (r, align) = get_multi_align(df, drop_threshold, column_threshold)
    if align:
        file_name = os.path.join(output_dir, name+ ".MULTI.cc.phy")
        with open(file_name,"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(align)
    return r

def get_morphological_df(name):
    df = pd.read_csv(os.path.join(input_dir, name + "/cldf/values.csv"), low_memory=False)
    columns = ["Language_ID", "Parameter_ID", "Value"]
    new_df = df
    for column in df.columns:
        if not column in columns:
            new_df = new_df.drop(column, axis=1)
    df = new_df

    indefinite_values = ["-", "?", "N/A", "X", "None"]
    for value in indefinite_values:
        df.replace(value, np.NaN)
    df = df[df.Value == df.Value]
    df = df.rename(columns={'Parameter_ID': 'Char_ID'})
    return df

def convert_morphological(name):
    df = get_morphological_df(name)
    align = get_bin_align(df)
    file_name = os.path.join(output_dir, name+ ".phy")
    with open(file_name,"w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(align)


def convert_morphological_multivalue(name, drop_threshold, column_threshold):
    df = get_morphological_df(name)
    (r, align) = get_multi_align(df, drop_threshold, column_threshold)
    if align:
        file_name = os.path.join(output_dir, name+ ".MULTI.ms.phy")
        with open(file_name,"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(align)
    return r





def convert_cldf_cognate_multi():
    column_threshold = 0
    drop_threshold = 0.1
    res_file = open("temp/cldf_cognate_multi.meta", "w+")
    res_file.write("name,num_languages,num_chars,multi_cells,multi_chars,high_multi_chars,high_state_chars,dropped_chars,max_states,converted\n")
    with os.scandir(input_dir) as it:
        for entry in it:
            if not entry.is_dir():
                continue
            p = os.path.join(input_dir, entry.name + "/cldf/")
            res = None
            if os.path.isfile(p + "forms.csv") and os.path.isfile(p + "cognates.csv"):
                print(entry.name)
                res = convert_cognate_multivalue(entry.name, drop_threshold, column_threshold)
            if res:
                string = (",").join([str(item) for key,item in res.items()])
                res_file.write(entry.name + "," + string + "\n")

def convert_cldf_morpho_multi():
    column_threshold = 0
    drop_threshold = 0.1
    res_file = open("temp/cldf_morpho_multi.meta", "w+")
    res_file.write("name,num_languages,num_chars,multi_cells,multi_chars,high_multi_chars,high_state_chars,dropped_chars,max_states,converted\n")
    with os.scandir(input_dir) as it:
        for entry in it:
            if not entry.is_dir():
                continue
            p = os.path.join(input_dir, entry.name + "/cldf/")
            res = None
            if os.path.isfile(p + "values.csv"):
                print(entry.name)
                res = convert_morphological_multivalue(entry.name, drop_threshold, column_threshold)
            if res:
                string = (",").join([str(item) for key,item in res.items()])
                res_file.write(entry.name + "," + string + "\n")




def convert_cldf_bin():
    with os.scandir(input_dir) as it:
        for entry in it:
            if not entry.is_dir():
                continue
            p = os.path.join(input_dir, entry.name + "/cldf/")
            if os.path.isfile(p + "forms.csv") and os.path.isfile(p + "cognates.csv"):
                print(entry.name)
                convert_cognate(entry.name)
            elif os.path.isfile(p + "values.csv"):
                print(entry.name)
                convert_morphological(entry.name)



def build_language_families():
    d = "raw/OSF/language_families"
    outfile = open("temp/language.families", "w+")
    with os.scandir(d) as it:
        for entry in it:
            if not entry.is_file() or not entry.name.endswith("cc.phy"):
                continue
            family = entry.name.split('.')[0]
            align = AlignIO.read(os.path.join(d, entry.name), "phylip-relaxed")
            languages = [rec.id.split('.')[-1] for rec in align]
            #languages = [rec.id for rec in align]
            outfile.write(family + ":")
            outfile.write(','.join(languages))
            outfile.write("\n")


def read_language_families():
    lines = open("temp/language.families", "r").read().split("\n")[:-1]
    family_dict = {}
    for line in lines:
        data = line.split(':')
        family_dict[data[0]] = data[1].split(',')
    return family_dict



def convert_osf_multi():
    column_threshold = 0.1
    drop_threshold = 0.1
    family_dict = read_language_families()
    df = pd.read_csv("raw/OSF/asjp17Clustered.csv")
    groups = df.groupby(df["glot_fam"])
    #family_matching(df, family_dict)
    res_file = open("temp/osf_multi.meta", "w+")
    res_file.write("name,num_languages,num_chars,multi_cells,multi_chars,high_multi_chars,high_state_chars,dropped_chars,max_states,converted\n")
    for family, sub_df in groups:
        if not family in family_dict:
            continue #too small

        sub_df = sub_df.rename(columns={'doculect': 'Language_ID', 'concept' : 'Char_ID', 'cClass' : 'Value'})
        (r, align) = get_multi_align(sub_df, drop_threshold, column_threshold)
        if align:
            file_name = os.path.join(output_dir, family+ ".MULTI.cc.phy")
            with open(file_name,"w+") as f:
                writer = RelaxedPhylipWriter(f)
                writer.write_alignment(align)
        string = (",").join([str(item) for key,item in r.items()])
        res_file.write(family + "," + string + "\n")


def convert_cldf_bin_catg(path, name):
    df = get_df_gerhard_csv(path)
    outfile = open("catg/" + name + ".catg", "w+")
    write_bin_catg(df, outfile)



#convert_cldf_cognate_multi()
#convert_cldf_morpho_multi()
#convert_osf_multi()

convert_cldf_bin_catg("raw/dunnielex_raw.csv", "dunnielex")





