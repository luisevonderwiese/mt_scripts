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

def aic(name, model):
    if model.startswith("BIN"):
        alignment_name = "morph_alignments/bin/" + name + ".BIN.phy"
    else:
        alignment_name = "morph_alignments/multi/" + name + ".phy"
    os.system('./../../tools/raxml-ng/build/bin/raxml-ng --search1 --msa ' + alignment_name +
            ' --threads 2 --model ' + model + ' --prefix foo --nofiles' +
              " --opt-branches off >out.txt")
    l_file = open('out.txt', 'r')
    lines = l_file.readlines()
    aic = []
    for line in lines:
        if(line.startswith('AIC score:')):
            aic.append(line.split(" ")[2].strip())
            aic.append(line.split(" ")[6].strip())
            aic.append(line.split(" ")[10].strip())
    os.remove("out.txt")
    aic_string = (',').join(aic)
    return aic_string


def calculate_aic_scores():
    outfile = open("temp/aic.scores", "w+")
    outfile.write("name|gtr|mk|bin\n")
    num_states_dict = read_num_states()
    with os.scandir("morph_alignments/multi") as it:
        for entry in it:
            if not entry.is_file():
                continue
            if not entry.name in num_states_dict:
                continue
            num_states = num_states_dict[entry.name]
            name = entry.name.split('.')[0]
            print(name)
            binmodel = "BIN"
            gtrmodel = "MULTI" + str(num_states) + "_GTR"
            mkmodel = "MULTI" + str(num_states) + "_MK"
            outfile.write(name)
            outfile.write('|')
            outfile.write(aic(name, gtrmodel))
            outfile.write('|')
            outfile.write(aic(name, mkmodel))
            outfile.write('|')
            outfile.write(aic(name, binmodel))
            outfile.write('\n')

calculate_aic_scores()
