from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os
import shutil
import re

states = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "!", "\"", "#", "$", "%", "&", "'", "(", ")", "*", "+",
",", "/", ":", ";", "<", "=", ">", "@", "[", "\\", "]", "^", "_", "{", "|", "}", "~"]
idx_for_state = {}
for i, state in enumerate(states):
    idx_for_state[state] = i


def convert_to_binary(align):
    c_count = len(align[0].seq)
    column_maxima = [0 for _ in range(c_count)]
    for record in align:
        sequence = record.seq
        for c in range(c_count):
            if(len(sequence) == c_count):
                value = sequence[c]
            else:
                value = sequence[c]['d'][0]
            if value == '-' or value == '?':
                continue
            column_maxima[c] = max(idx_for_state[value], column_maxima[c])
    new_sequences = []
    for record in align:
        new_sequence = ""
        sequence = record.seq
        for c in range(c_count):
            block_length = column_maxima[c] + 1
            if(len(sequence) == c_count):
                value = sequence[c]
            else:
                value = sequence[c]['d'][0]
            if value == '-' or value == '?':
                new_value = "".join([value for _ in range(block_length)])
            else:
                chars = ['0' for _ in range(block_length)]
                chars[idx_for_state[value]] = '1'
                new_value = "".join(chars)
            new_sequence += new_value
        new_sequences.append(new_sequence)
    new_records = [SeqRecord(new_sequences[i], id=align[i].id) for i in range(len(align))]
    return MultipleSeqAlignment(new_records, annotations={}, column_annotations={})

def convert_multivalue_phylip(input_dir, output_dir, name):
    align = AlignIO.read(os.path.join(input_dir, name + ".phy"), "phylip-relaxed")
    check_alignment_for_special_chars(os.path.join(input_dir, name + ".phy"))
    align = convert_to_binary(align)
    with open(os.path.join(output_dir, name + ".BIN.phy"),"w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(align)


def convert_mulitvalue_csv(input_dir, output_dir, name):
    indefinite_values = ["-", "?", "N/A", "X", "None"]
    with open(os.path.join(input_dir, name + ".csv"), "r") as values_file:
        lines = values_file.readlines()
    values_for_params = {}
    data_len = len(lines[0].split(","))
    for i in range(1, len(lines)):
        data = lines[i].split(",")
        if len(data) < data_len:
            continue
        param_id = data[2]
        value = data[3]
        if value in indefinite_values:
            continue
        if not param_id in values_for_params.keys():
            values_for_params[param_id] = set()
        values_for_params[param_id].add(value)
    translate = {}
    max_states = 0
    for (param_id, values) in values_for_params.items():
        values = list(values)
        values.sort()
        #print(values)
        max_states = max(len(values), max_states)
        param_translate = {}
        for (i, value) in enumerate(values):
            param_translate[value] = states[i]
        translate[param_id] = param_translate
    print(name + " max_states " + str(max_states))
    matrix = {}
    for i in range(1, len(lines)):
        data = lines[i].split(",")
        if len(data) < data_len:
            continue
        language_id = data[1]
        param_id = data[2]
        if not language_id in matrix:
            matrix[language_id] = {}
        if param_id in matrix[language_id]:
            print(param_id)
            print(language_id)
            assert(not param_id in matrix[language_id])
        value = data[3]
        if value in indefinite_values:
            matrix[language_id][param_id] = "-"
        else:
            matrix[language_id][param_id] = translate[param_id][value]
    all_param_ids = set()
    for language_id in matrix.keys():
        all_param_ids = all_param_ids.union(set(matrix[language_id].keys()))
    for language_id in matrix.keys():
        for param_id in all_param_ids:
            if not param_id in matrix[language_id].keys():
                matrix[language_id][param_id] = "-"
    all_param_ids = list(all_param_ids)
    all_param_ids.sort()
    sequences = {}
    for language_id in matrix.keys():
        s = ""
        for param_id in all_param_ids:
            s += matrix[language_id][param_id]
        sequences[language_id] = s
    records = [SeqRecord(sequences[language_id], id=language_id) for language_id in matrix.keys()]
    align = MultipleSeqAlignment(records, annotations={}, column_annotations={})
    file_name = os.path.join(output_dir, name+ ".phy")
    with open(file_name,"w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(align)
    #shutil.copyfile(file_name, os.path.join("conversion_inputs/multivalue_phylip/", name +  ".phy"))
