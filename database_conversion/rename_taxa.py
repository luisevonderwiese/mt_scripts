import os
import pandas as pd

from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def rename_taxa():
    dir = "taxa_name_mappings"
    with os.scandir(dir) as it:
        for entry in it:
            if not entry.is_file():
                continue


            lines = open(os.path.join(dir, entry.name), "r").readlines()[1:]
            d={}
            for line in lines:
                data = line.split(",")
                if (len(data) != 2):
                    print(line)
                    continue
                d[data[0]] = data[1]

            alignment_name = ".".join(entry.name.split('.')[:-1]) + ".phy"
            alignment_path = os.path.join("alignment_downloads", alignment_name)
            align = AlignIO.read(alignment_path, "phylip-relaxed")
            new_records = []
            for rec in align:
                if rec.id in d:
                    new_records.append(SeqRecord(rec.seq, id=d[rec.id]))
                else:
                    new_records.append(SeqRecord(rec.seq, rec.id))
            new_align = MultipleSeqAlignment(new_records, annotations={}, column_annotations={})
            new_alignemnt_path = os.path.join("alignment_uploads/alignments", alignment_name)
            with open(new_alignemnt_path,"w+") as f:
                writer = RelaxedPhylipWriter(f)
                writer.write_alignment(new_align)


rename_taxa()
