from Bio import AlignIO
from ete3 import Tree
import os

alignment_dir = "duplicate_check/"
mapping = {"NIIGATA": "NIGATA",
           "IBARAKI": "IBARAGI",
           "HACHIJYO": "HACHIJO",
           "NORTHALTAY" : "NORTHALTAI",
           "SOUTHALTAY" : "SOUTHALTAI",
           "BARABATATAR" : "BARABA",
           "CUMAN" : "CODEXCUMANICUS",
           "UYGHUR" : "UIGHUR"
           }



def taxa_name_unification_tree(tree):
    for leaf in tree.iter_leaves():
        taxa_id = leaf.name
        taxa_id = taxa_id.upper()
        taxa_id = taxa_id.replace('_', '')
        if taxa_id in mapping.keys():
            taxa_id = mapping[taxa_id]
        leaf.name = taxa_id
    return tree

def taxa_name_unification(align):
    for record in align:
        taxa_id = record.id
        taxa_id = taxa_id.upper()
        taxa_id = taxa_id.replace('_', '')
        if taxa_id in mapping.keys():
            taxa_id = mapping[taxa_id]
        record.id = taxa_id
    return align

def replace_indef_zero(align):
    for rec in align:
        rec.seq = rec.seq.replace('-', '0')
    return align

def check_duplicates(alignment_name1, alignment_name2):
    align1 = AlignIO.read(alignment_dir + alignment_name1, "phylip-relaxed")
    align2 = AlignIO.read(alignment_dir + alignment_name2, "phylip-relaxed")
    #align1 = replace_indef_zero(align1)
    num_sites1 = align1.get_alignment_length()
    num_sites2 = align2.get_alignment_length()
    if num_sites1 != num_sites2:
        print(num_sites1)
        print(num_sites2)
        print("!!! Different number of sites")
        return False
    #align1 = taxa_name_unification(align1)
    #align2 = taxa_name_unification(align2)
    taxa1 = set([align1[i].id for i in range(len(align1))])
    taxa2 = set([align2[i].id for i in range(len(align2))])
    common_taxa = taxa1.intersection(taxa2)
    taxa1_only = taxa1.difference(common_taxa)
    taxa2_only = taxa2.difference(common_taxa)
    if (len(taxa1_only) != 0 or len(taxa2_only) != 0):
        print(taxa1_only)
        print(taxa2_only)
        print("!!! Different Taxa")
        return False
    align1.sort(key=lambda x: x.id)
    align2.sort(key=lambda x: x.id)
    columns1 = [align1[:, i] for i in range(num_sites1)]
    columns2 = [align2[:, i] for i in range(num_sites2)]

    L1 = [ (columns1[i],i) for i in range(len(columns1)) ]
    L1.sort()
    columns1, permutation1 = zip(*L1)

    L2 = [ (columns2[i],i) for i in range(len(columns2)) ]
    L2.sort()
    columns2, permutation2 = zip(*L2)

    j = 0
    i = 0
    while i < len(columns1) and j < len(columns2):
        if columns1[i] == columns2[j]:
            i+=1
            j+=1
            continue
        if columns1[i] < columns2[j]:
            print("Unmatch 1")
            print(columns1[i])
            #print(permutation1[i])
            i+=1
            return False
        if columns1[i] > columns2[j]:
            print("Unmatch 2")
            print(columns2[j])
            #print(permutation2[j])
            j+=1
            return False
    return True


def rf_distance(alignment_name1, alignment_name2):
    print(alignment_name1)
    print(alignment_name2)
    t1 = Tree(alignment_dir + alignment_name1 + ".raxml.bestTree")
    t2 = Tree(alignment_dir + alignment_name2 + ".raxml.bestTree")
    t1 = taxa_name_unification_tree(t1)
    t2 = taxa_name_unification_tree(t2)
    rf, max_rf, common_leaves, parts_t1, parts_t2,discard_t1, discard_t2 = t1.robinson_foulds(t2, unrooted_trees = True)
    print(rf)
    print(max_rf)
    print(common_leaves)
    #print(parts_t1)
    #äprint(parts_t2)
    #print(discard_t1)
    #print(discard_t2)
    if max_rf == 0:
        print("?!")
        return 0
    rel_rf = rf/max_rf
    print("RF Distance of ML Trees: " + str(rel_rf))
    return rel_rf




#Grollemund_et_al ist subalign
#Bei den anderen stimmen zumindest die Taxa überein, (birchall un dlee2015 musste ich die von phlorest umbenennen)
#Gleich sind aber keine

pnd = {}
pnd["dyenindoeuropean"] = ["gray_and_atkinson2003"]
pnd["grollemundbantu"] = ["grollemund_et_al2015", "koile_et_al2022"]
pnd["kitchensemitic"] = ["kitchen_et_al2009"]
pnd["powerma"] = ["power_et_al2020"]
#pnd["iecor"] = "heggarty_et_al_subm" (iecor not in lexibank)
pnd["birchallchapacuran"] = ["birchall_et_al2016"]
pnd["gerarditupi"] = ["gerardi_and_reichert2021"]
pnd["utoaztecan"] = ["greenhill_et_al_subm"]
pnd["dravlex"] = ["kolipakam_et_al2018"]
pnd["leekoreanic"] = ["lee2015"]
pnd["leejaponic"] = ["lee_and_hasegawa2011"]
pnd["leeainu"] = ["lee_and_hasegawa2013"]
pnd["nagarajakhasian"] = ["nagaraja_et_al2013"]
pnd["robinsonap"] = ["robinson_and_holton2012"]
pnd["sagartst"] = ["sagart_et_al2019"]
pnd["savelyevturkic"] = ["savelyev_and_robbeets2020"]
#pnd["walkerarawakan"] = "walker_and_ribeiro2011" (no cognates.csv in lexibank)
pnd["mcelhanonhuon"] = ["greenhill2015"]

for (lexiname, phlorestnames) in pnd.items():
    lexiname="selfmade/" + lexiname + ".phy"
    for phlorestname in phlorestnames:
        print("\n" + phlorestname)
        phlorestname = "phlorest/" + phlorestname + ".cc.phy"
        check_duplicates(lexiname, phlorestname)
    if len(phlorestnames) > 1:
        phlorestname0 = "phlorest/" + phlorestnames[0] + ".cc.phy"
        phlorestname1 = "phlorest/" + phlorestnames[1] + ".cc.phy"
        check_duplicates(phlorestname0, phlorestname1)


