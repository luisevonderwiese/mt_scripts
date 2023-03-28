from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os
import shutil
import re



def convert_nexus_to_phylip(input_dir, output_dir, alignment_name):
    align = AlignIO.read(input_dir + alignment_name + ".nex", "nexus")
    outfile_name = os.path.join(output_dir, alignment_name +  ".phy")
    with open(outfile_name,"w+") as f:
        writer = RelaxedPhylipWriter(f)
        writer.write_alignment(align)

def proper_nex_heads(input_dir, output_dir, alignment_name):
    def get_head(n_taxa, n_char):
        return "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=" + str(n_taxa) + " NCHAR=" + str(n_char) + ";\n\tFORMAT DATATYPE=Standard SYMBOLS= \"01\" MISSING=? GAP= -;\nMATRIX\n"
    end = "\n ;\nend;\n"

    content = open(os.path.join(input_dir, alignment_name), "r").read()
    content = "\n".join(content.split("\n")[:-2]) + end
    parts = content.split("MATRIX")
    line = parts[0].split("\n")[2].split(" ")
    n_char = line[4][:-1]
    n_taxa = line[1].split('=')[1]
    align = parts[1][2:]
    new_content = get_head(n_taxa, n_char) + align
    open(os.path.join(output_dir, alignment_name), "w+").write(new_content)


def convert():
    output_dir = "conversion_inputs/nexus/"
    input_dir = "conversion_inputs/nexus_wrong/"
    with os.scandir(input_dir) as it:
        for entry in it:
            if not entry.is_file():
                continue
            proper_nex_heads(input_dir, output_dir, entry.name)
    input_dir = "conversion_inputs/nexus/"
    output_dir = "alignment_uploads/alignments/"
    with os.scandir(input_dir) as it:
        for entry in it:
            if not entry.is_file():
                continue
            print(entry.name)
            convert_nexus_to_phylip(input_dir, output_dir, ".".join(entry.name.split(".")[:-1]))
    with os.scandir(output_dir) as it:
        for entry in it:
            if not entry.is_file():
                continue
            try:
                align = AlignIO.read(os.path.join(output_dir, entry.name), "phylip-relaxed")
            except:
                print(entry.name + " padding required!")
                continue
            for record in align:
                if not re.match('^[a-zA-Z0-9_\\-\\.]+$', record.id):
                    print(entry.name + ": illegal char in taxa name " + record.id)


convert()

