import pathlib
from Bio import SeqIO, AlignIO
import re

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


AMBIG_CODE = {
    # Nucleotides
    "{CT}": "Y",
    "{TC}": "Y",
    "{AG}": "R",
    "{GA}": "R",
    "{AT}": "W",
    "{TA}": "W",
    "{CG}": "S",
    "{GC}": "S",
    "{TG}": "K",
    "{GT}": "K",
    "{CA}": "M",
    "{AC}": "M",

    "{AGT}": "D",
    "{ATG}": "D",
    "{GTA}": "D",
    "{GAT}": "D",
    "{TAG}": "D",
    "{TGA}": "D",

    "{ACG}": "V",
    "{AGC}": "V",
    "{CGA}": "V",
    "{CAG}": "V",
    "{GCA}": "V",
    "{GAC}": "V",

    "{ACT}": "H",
    "{ATC}": "H",
    "{CTA}": "H",
    "{CAT}": "H",
    "{TAC}": "H",
    "{TCA}": "H",

    "{CGT}": "B",
    "{CTG}": "B",
    "{GTC}": "B",
    "{GCT}": "B",
    "{TCG}": "B",
    "{TGC}": "B",

    "{ACGT}": "N",
    "{ACTG}": "N",
    "{AGCT}": "N",
    "{AGTC}": "N",
    "{ATCG}": "N",
    "{ATGC}": "N",
    "{CAGT}": "N",
    "{CATG}": "N",
    "{CGAT}": "N",
    "{CGTA}": "N",
    "{CTAG}": "N",
    "{CTGA}": "N",
    "{GACT}": "N",
    "{GATC}": "N",
    "{GCAT}": "N",
    "{GCTA}": "N",
    "{GTAC}": "N",
    "{GTCA}": "N",
    "{TACG}": "N",
    "{TAGC}": "N",
    "{TCAG}": "N",
    "{TCGA}": "N",
    "{TGAC}": "N",
    "{TGCA}": "N",

    # Amino Acids
    "{DN}": "B",
    "{ND}": "B",
    "{EQ}": "Z",
    "{QE}": "Z",

    # other weird characters representing a gap
}

GAP_MAP = {
    ":": "-",
    "#": "-",
    "?": "-",
    "~": "-",
    ".": "-"
}


def get_data_type(nexus_metadata):
    data_type = ""

    for line in nexus_metadata:
        line = line.strip().lower()
        if "datatype" in line:
            # FORMAT DATATYPE=DNA SYMBOLS= "A C G T" MISSING=? GAP= -;
            if "dna" in line or "nucleotide" in line or "rna" in line:
                data_type = "dna"
            elif "protein" in line:
                data_type = "aa"
            elif "standard" in line:
                data_type = "standard"
            else:
                # in this case some chars are numbers
                data_type = "other"

    return data_type


def get_matrix_sections_and_metadata(nexus_file):
    content = open(nexus_file).readlines()
    matrices = []
    metadata = []
    partitions = []

    seen_matrix_begin = False
    seen_matrix_end = False

    seen_part_begin = False

    matrix_begin = None
    metadata_begin = None
    partition_begin = None

    matrix_counter = 0

    for i, l in enumerate(content):
        l = l.strip()
        l = l.upper()

        if l.startswith("MATRIX"):
            seen_matrix_begin = True
            seen_matrix_end = False
            matrix_begin = i
            # begin of the matrix means the metadata information section ended in the line above
            metadata.append((metadata_begin, i - 1))
        elif seen_matrix_begin and l.startswith("END;"):
            matrices.append((matrix_begin, i))
            # for each matrix we find, we add a dummy partition
            # in case we find a partition definition later, we replace this item
            partitions.append((None, None))
            seen_matrix_begin = False
            matrix_counter += 1
            matrix_begin = None
            seen_matrix_end = True
        elif l.startswith("BEGIN CHARACTERS"):
            metadata_begin = i
        elif l.startswith("BEGIN SETS") and seen_matrix_end:
            seen_part_begin = True
            seen_matrix_end = False
            partition_begin = i
        elif seen_part_begin and l.startswith("END;"):
            seen_part_begin = False
            # replace the dummy partition for the matrix with the correct partition section
            try:
                partitions[matrix_counter - 1] = (partition_begin, i)
            except:
                print(nexus_file, matrix_counter)
                assert 0

    assert len(matrices) == len(metadata) == len(partitions), f"{nexus_file}: {len(matrices)}, {len(metadata)}, {len(partitions)}"

    return matrices, metadata, partitions


def format_taxon_name(taxon):
    taxon = re.sub(f"\s+", " ", taxon)
    taxon = taxon.strip().strip("'")
    taxon = taxon.replace(" ", "_")
    taxon = taxon.replace(";", "")
    taxon = taxon.replace(":", "")
    taxon = taxon.replace(",", "")
    taxon = taxon.replace("(", "_")
    taxon = taxon.replace(")", "_")
    taxon = taxon.replace("'", "")
    return taxon


def get_model(msa_file, data_type):
    if data_type == "dna":
        return "GTR+G"
    elif data_type == "aa":
        return "LG+G"
    elif data_type == "standard":
        msa = AlignIO.read(msa_file, format="phylip-relaxed")
        unique_states = set()
        for seq in msa:
            seq = seq.seq.replace("-", "")
            unique_states = unique_states.union(set(seq))
        # the number of unique states is irrelevant for RAxML-NG, it only cares about the max state value...
        num_states = int(max(unique_states)) + 1
        return f"MULTI{num_states}_GTR"
    else:
        raise ValueError("Unsupported data type: ", data_type)


def convert_nexus_to_phylip(nexus_matrix):
    """
    Converts the MSA section in Nexus format to a MSA in phylip-relaxed format
    :param nexus_matrix: List of raw lines out of the nexus file
    :return: List of SeqRecord objects containing the rows of the MSA
    """
    sequences = []
    taxons = []
    for m in nexus_matrix:
        m = m.strip()
        if m.startswith("["):
            continue
        try:
            taxon, seq = m.rsplit(" ", 1)
        except:
            continue

        taxon = format_taxon_name(taxon)

        # hack to prevent duplicate taxon names
        ct = 1
        while taxon in taxons:
            taxon += str(ct)
            ct += 1
        taxons.append(taxon)
        seq = seq.strip()

        for key, value in AMBIG_CODE.items():
            seq = seq.replace(key, value)
        for key, value in GAP_MAP.items():
            seq = seq.replace(key, value)
        sequences.append(SeqRecord(Seq(seq), id=taxon, name=taxon))

    return sequences


def convert_nexus_morph_data_to_phylip(nexus_matrix):
    """
    Converts the MSA section in Nexus format to a MSA in phylip-relaxed format
    This method is specifically for morphological data
    :param nexus_matrix: List of raw lines out of the nexus file
    :return: List of SeqRecord objects containing the rows of the MSA
    """
    sequences = []
    taxons = []
    for m in nexus_matrix:
        m = m.strip()
        if m.startswith("["):
            continue
        try:
            taxon, seq = m.rsplit(" ", 1)
        except:
            continue

        taxon = format_taxon_name(taxon)

        # hack to prevent duplicate taxon names
        ct = 1
        while taxon in taxons:
            taxon += str(ct)
            ct += 1
        taxons.append(taxon)
        seq = seq.strip()

        # replace all occurences of ambiguity, e.g. {01} with a gap char
        seq = re.sub(r"{\d+}", "-", seq)

        for key, value in GAP_MAP.items():
            seq = seq.replace(key, value)

        morph_regex = re.compile(r"([0-9\-\{\}]+)")
        candidates = re.findall(morph_regex, seq)

        morph_seqs = []
        for candidate_seq in candidates:
            if set(candidate_seq) != {"-"}:
                morph_seqs.append(candidate_seq)

        if len(morph_seqs) != 1:
            raise ValueError(f"Found multiple or no possible sequences: {morph_seqs}")

        seq = morph_seqs[0]

        sequences.append(SeqRecord(Seq(seq), id=taxon, name=taxon))

    return sequences


def get_partitions(partition_info):
    partitions = []
    reg = re.compile(r"charset\s+([^\s]+).*")

    for line in partition_info:
        line = line.strip().lower()
        if not line.startswith("charset"):
            continue
        # charset wg (CHARACTERS = 'Heliconiini, individuals, full dataset') =  4035-4408;
        m = re.match(reg, line)
        name = ""
        if m:
            name = m.groups()[0]

        parts = line.rsplit("=", 1)[1].strip(";").strip().split(" ")

        partitions.append((name, parts))

    return partitions


def convert_and_save_nexus(nexus_file):
    filename_base = pathlib.Path(nexus_file).name

    content = open(nexus_file).readlines()
    matrices, metadata, partitions = get_matrix_sections_and_metadata(nexus_file)

    for i, (
        (matrix_begin, matrix_end),
        (metadata_begin, metadata_end),
        (partition_begin, partition_end),
    ) in enumerate(zip(matrices, metadata, partitions)):
        filename = filename_base + f"_{i}.phy"

        metadata = content[metadata_begin:metadata_end]
        data_type = get_data_type(metadata)

        matrix = content[matrix_begin:matrix_end]
        if data_type in ["other", "standard"]:
            try:
                sequences = convert_nexus_morph_data_to_phylip(matrix)
            except ValueError as e:
                print("Failed to convert nexus morph data to phylip ", nexus_file, i, "error is: ", e)
                return

        else:
            sequences = convert_nexus_to_phylip(matrix)

        if len(sequences) == 0:
            print(f"Found no sequences in matrix {i} of file {nexus_file}")
            continue

        n_sites = len(sequences[0].seq)

        if n_sites <= 1:
            print("Sequences contain only 1 character, ", nexus_file, i)
            continue

        all_n_sites = set([len(seq.seq) for seq in sequences])
        if len(all_n_sites) > 1:
            print("Sequences are of different lengths: ", all_n_sites, nexus_file, i)
            continue

        try:
            SeqIO.write(sequences, filename, "phylip-relaxed")
        except Exception as e:
            print("Error writing phylip ", nexus_file, i, "error", e)
            continue

        try:
            model = get_model(filename, data_type)
        except ValueError as e:
            print("Error reading phylip ", filename)
            print(e)
            continue

        partition_filename = filename_base + f"_{i}.part"

        partition_info = content[partition_begin:partition_end]
        partitions = get_partitions(partition_info)

        part_string = ""

        for (name, parts) in partitions:
            part_string += f"{model}, {name} = {', '.join(parts)}\n"

        open(partition_filename, "w").write(part_string)


convert_and_save_nexus("example.nex")