import re
import json
from gfagraphs import Graph
from itertools import pairwise
from collections import defaultdict
from pathlib import Path
from datetime import datetime


def parse_snarls(json_path: str):
    # Load snarl json file and return a list of boundaries and a list of parents with the same order
    # The index of a element in a list corresponds to the index of the structure
    # The boundaries and parents are stored as tuples of (start, end) node ids, with start < end
    # If there is no parent, the parent is set to None
    boundaries = []
    parents = []

    with open(json_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            data = json.loads(line)

            start = int(data["start"]["node_id"])
            end = int(data["end"]["node_id"])
            if start > end:
                temp = start
                start = end
                end = temp
            key = (str(start), str(end))

            parent = None
            if "parent" in data:
                p = data["parent"]
                p_start = int(p["start"]["node_id"])
                p_end = int(p["end"]["node_id"])
                if p_start > p_end:
                    temp = p_start
                    p_start = p_end
                    p_end = temp
                parent = (str(p_start), str(p_end))

            boundaries.append(key)
            parents.append(parent)

    return boundaries, parents


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

def write_report(
    output_dir,
    input_gfa,
    max_len,
    min_variant,
    save_subgraphs,
    compression,
    vg_snarls_timestamp,
    bubbles,
    start_time,
    elapsed,
    peak_ram_gb,
):

    report = Path(output_dir) / "simplification_report.txt"

    initial = compression["initial_nodes"]
    final = bubbles["after_small_variants"]

    total_removed = initial - final

    # Detailed execution times
    graph_load_time = compression["graph_load_time"] - start_time
    linear_compression_time = compression["linear_compression_time"] - compression["graph_load_time"]
    snp_mnp_compression_time = compression["snp_mnp_compression_time"] - compression["linear_compression_time"]
    temporary_saving_time = compression["temporary_saving_time"] -  compression["snp_mnp_compression_time"]
    vg_snarls_time = vg_snarls_timestamp - compression["temporary_saving_time"]
    bubbles_load_and_index_time = bubbles["bubbles_load_and_index_time"] - vg_snarls_timestamp
    compressed_graph_load_time = bubbles["graph_load_time"] - bubbles["bubbles_load_and_index_time"]
    node_and_size_extraction_time = bubbles["node_and_size_extraction_time"] - bubbles["graph_load_time"]
    bubbles_selection_time = bubbles["bubbles_selection_time"] - bubbles["node_and_size_extraction_time"]
    bubbles_abstraction_time = bubbles["bubbles_abstraction_time"] - bubbles["bubbles_selection_time"]
    if bubbles["subgraph_build_time"] != 0:
        subgraph_build_time = bubbles["subgraph_build_time"] - bubbles["bubbles_abstraction_time"]
        subgraph_save_time = bubbles["subgraph_save_time"] - bubbles["subgraph_build_time"]
        main_graph_build_time = bubbles["main_graph_build_time"] - bubbles["subgraph_save_time"]
    else:
        subgraph_build_time = 0
        subgraph_save_time = 0
        main_graph_build_time = bubbles["main_graph_build_time"] - bubbles["bubbles_abstraction_time"]
    graph_save_time = bubbles["graph_save_time"] - bubbles["main_graph_build_time"]



    with open(report, "w") as f:

        f.write("=" * 60 + "\n")
        f.write("PGSIMPLIFY REPORT\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Date                 : {datetime.now()}\n")
        f.write(f"Input graph          : {input_gfa}\n")
        f.write(f"Output directory     : {output_dir}\n\n")

        f.write("Parameters\n")
        f.write("-" * 30 + "\n")
        f.write(f"Maximum MNP length   : {max_len}\n")
        f.write(f"Minimum variant size : {min_variant}\n")
        f.write(f"Save subgraphs       : {save_subgraphs}\n\n")

        f.write("Compression summary\n")
        f.write("-" * 30 + "\n")

        f.write(f"Initial node number                  : {compression['initial_nodes']}\n\n")

        f.write("Non-branching path compression\n")
        f.write(f"    Removed nodes                    : {compression['removed_non_branching']}\n")
        f.write(f"    Remaining nodes                  : {compression['after_non_branching']}\n\n")

        f.write("SNP/MNP compression\n")
        f.write(f"    Removed nodes                    : {compression['removed_snp']}\n")
        f.write(f"    Remaining nodes                  : {compression['after_snp']}\n\n")

        f.write("Small variant simplification\n")
        f.write(f"    Bubbles detected                  : {bubbles['nb_bubbles_detected']}\n")
        f.write(f"    Bubbles simplified                : {bubbles['nb_bubbles_simplified']}\n")
        f.write(f"    Removed nodes                    : {bubbles['removed_small_variants']}\n")
        f.write(f"    Remaining nodes                  : {bubbles['after_small_variants']}\n\n")

        f.write("Global statistics\n")
        f.write("-" * 30 + "\n")
        f.write(f"Initial nodes                        : {initial}\n")
        f.write(f"Final nodes                          : {final}\n")
        f.write(f"Total removed                        : {total_removed}\n")
        f.write(f"Global reduction                     : {100 * total_removed / initial:.2f} %\n\n")

        f.write("Execution\n")
        f.write("-" * 30 + "\n")
        f.write(f"Wall time                            : {elapsed:.2f} s\n")
        f.write(f"Peak RAM                             : {peak_ram_gb:.2f} GB\n")

        f.write("Detailed execution times\n")
        f.write("-" * 30 + "\n")
        f.write(f"Graph load time (gfagraphs)           : {graph_load_time:.2f} s\n")
        f.write(f"Linear compression time               : {linear_compression_time:.2f} s\n")
        f.write(f"SNP/MNP compression time              : {snp_mnp_compression_time:.2f} s\n")
        f.write(f"Temporary saving time (gfagraphs)     : {temporary_saving_time:.2f} s\n")
        f.write(f"Snarls detection time (vg snarls)     : {vg_snarls_time:.2f} s\n")
        f.write(f"Bubbles load and index time           : {bubbles_load_and_index_time:.2f} s\n")
        f.write(f"Compressed graph load time (gfagraphs): {compressed_graph_load_time:.2f} s\n")
        f.write(f"Node and size extraction time         : {node_and_size_extraction_time:.2f} s\n")
        f.write(f"Bubbles selection time                : {bubbles_selection_time:.2f} s\n")
        f.write(f"Bubbles abstraction time              : {bubbles_abstraction_time:.2f} s\n")
        if subgraph_build_time != 0:
            f.write(f"Subgraph build time                   : {subgraph_build_time:.2f} s\n")
            f.write(f"Subgraph save time (gfagraphs)        : {subgraph_save_time:.2f} s\n")
        f.write(f"Main graph build time                 : {main_graph_build_time:.2f} s\n")
        f.write(f"Main graph save time (gfagraphs)      : {graph_save_time:.2f} s\n")
