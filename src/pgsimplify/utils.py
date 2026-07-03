import re
from gfagraphs import Graph
from itertools import pairwise
from collections import defaultdict
from pathlib import Path


def compute_edge_orientation(graph):
    """
    Compute orientation of the nodes at each side of an edge for each path passing by this edge

    Parameter
    ----------
    graph : Graph
        The pangenome graph containing all edges and information to compute edges orientations
    
    Returns
    -------
    dict(list[(str,str)])
        Dictionnary (edge_begin, edge_end) -> orientation list
    """
    edge_orient = defaultdict(set)

    for path in graph.paths.values():

        for (u, ou), (v, ov) in pairwise(path["path"]):

            edge_orient[(u, v)].add(
                (ou.value, ov.value)
            )

    return edge_orient

def invert_orient(orient):
    """
    Invert an edge orientation list

    Parameter
    ----------
    orient : dict(list[(str,str)])
        Orientation list

    Returns
    -------
    dict(list[(str,str)])
        The inverted orientation list
    """
    inverted = set()

    for o1, o2 in orient:
        inverted.add((
            '-' if o2 == '+' else '+',
            '-' if o1 == '+' else '+'
        ))

    return inverted

def revcomp(string: str) -> str:
    """
    Compute reverse complement of a sequence
    Supports IUPAC ambiguity codes

    Parameter
    ----------
    string : str
        The sequence to reverse complement
    
    Returns
    -------
    str
        the reverse complement of input string
    """
    compl = {
        'A': 'T','C': 'G','G': 'C','T': 'A','N': 'N',
        "R": "Y","Y": "R","S": "S","W": "W",
        "K": "M","M": "K","B": "V","V": "B",
        "D": "H","H": "D"
    }
    return ''.join(compl[s] for s in string[::-1])

def get_substitutor(letters):
    """
    Gives IUPAC ambiguity code for a list of sequences to compress in one

    Parameter
    ----------
    letters : List(str)
        The sequences to compress

    Returns
    -------
    str
        The IUPAC ambiguity code corresponding to input sequences
    """
    out = ''
    for i in range(len(letters[0])):
        col = {s[i] for s in letters}

        if col == {'A','G'}: out += 'R'
        elif col == {'C','T'}: out += 'Y'
        elif col == {'G','T'}: out += 'K'
        elif col == {'A','C'}: out += 'M'
        elif col == {'C','G'}: out += 'S'
        elif col == {'A','T'}: out += 'W'
        elif col == {'A','C','G'}: out += 'V'
        elif col == {'A','C','T'}: out += 'H'
        elif col == {'A','G','T'}: out += 'D'
        elif col == {'C','G','T'}: out += 'B'
        elif col == {'A','C','G','T'}: out += 'N'
        else:
            out += next(iter(col))

    return out

def load_offsets(offset_file):
    """
    Load the dictionnary path -> offset
    Ex: #haplotype_name	offset
        hap1	42
        hap2	17
    
    Parameter
    ----------
    offset_file : 
        Path to offset file

    Returns
    -------
    dict[int]
        Offsets dictionnary
    """
    offsets = {}

    with open(offset_file, "r") as f:
        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            path_name, offset = line.split("\t")
            offsets[path_name] = int(offset)

    return offsets


def save_offsets(offsets, output_file):
    """
    Saves the dictionnary path -> offset in tab-delimited file

    Parameter
    ----------
    offsets : 
        Offsets dictionnary 
    output_file : 
        Directory to save the file 
    """
    with open(output_file, "w") as f:
        f.write(f"#haplotype_name\toffset\n")

        for haplotype_name, offset in offsets.items():
            f.write(f"{haplotype_name}\t{offset}\n")



def load_subgraphs(subgraphs_dir):
    """
    Load all subgraphs (gfa format) and index them by name
    
    Parameter
    ----------
    subgraphs_dir : str
        Directory where all subgraphs are stored

    Returns
    -------
    dict[Graph]
        Subgraph dictionnary indexed by their name
    """
    subgraphs = {}

    for gfa_file in sorted(Path(subgraphs_dir).glob("*.gfa")):

        match = re.search(r"(sg\d+)", gfa_file.stem)

        if match is None:
            continue

        sg_name = match.group(1)

        subgraphs[sg_name] = Graph(str(gfa_file), with_sequence=True)

    return subgraphs
