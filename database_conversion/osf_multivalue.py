
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import pandas as pd
from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os


def multistates_for_family(family, df):
    languages = list(set(df['doculect']))
    print(str(len(languages)) + " languages")
    matrix = {}
    for language in languages:
        matrix[language] = ""
    concepts = list(set(df['concept']))
    concepts.sort()
    num_concepts = len(concepts)
    print(str(num_concepts) + " concepts")
    colwise_multistates = []
    multistates = []
    for (c, concept) in enumerate(concepts):
        col_max = 0
        sub_df = df[df["concept"] == concept]
        concept_languages = list(set(sub_df['doculect']))
        classes = list(set(sub_df['cClass']))
        if (len(classes) > len(states)):
            continue
            #return (num_concepts, num_concepts, max_states, False)
        counts = sub_df['cClass'].value_counts()
        classes = sorted(classes, key=counts.get, reverse = True)
        for language in languages:
            if language not in concept_languages:
                multistates.append(0)
                continue
            sub_sub_df = sub_df[(sub_df["doculect"] == language)]
            class_count = len(sub_sub_df)
            col_max = max(col_max, class_count)
            multistates.append(class_count)
        colwise_multistates.append(col_max)
    return (multistates, colwise_multistates)


def multistate_statistic():
    family_dict = read_language_families()
    df = pd.read_csv("raw/OSF/asjp17Clustered.csv")
    groups = df.groupby(df["glot_fam"])
    multistates = []
    multistates_colwise = []
    for family, sub_df in groups:
        if not family in family_dict:
            continue #too small
        res = multistates_for_family(family, sub_df)
        multistates += res[0]
        multistates_colwise += res[1]
    plt.xlabel("Number of cognate classes")
    plt.ylabel("Number of Elements in (Languages x Concepts)")
    plt.hist(multistates)
    plt.savefig("temp/osf_multistate_statistics.png")
    plt.show()

    plt.xlabel("Maximum number of cogante classes per language")
    plt.ylabel("Number of Concepts")
    plt.hist(multistates_colwise, bins=20)
    plt.savefig("temp/osf_multistate_statistics_colwise.png")



#convert()
multistate_statistic()


