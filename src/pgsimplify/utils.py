import re
import json
import subprocess
from gfagraphs import Graph
from itertools import pairwise
from collections import defaultdict
from pathlib import Path
from datetime import datetime


def compute_snarls(tmpdir: Path):
    """
    Calls vg snarls and outputs snarls in json format
    
    Parameter
    ----------
    tmpdir : Path
        path to compressed graph in which we want to detect snarls
    """
    gfa_file = tmpdir / "compressed_graph.gfa"
    vg_file = tmpdir / "graph.vg"
    snarls_file = tmpdir / "graph.snarls"
    json_file = tmpdir / "bubbles.json"

    # GFA graph file -> VG graph file
    with open(vg_file, "wb") as out:
        subprocess.run(
            ["vg", "convert", "-f", "-W", str(gfa_file)],
            stdout=out,
            check=True,
        )

    # VG graph file -> binary snarls file
    with open(snarls_file, "wb") as out:
        subprocess.run(
            ["vg", "snarls", str(vg_file)],
            stdout=out,
            check=True,
        )

    # binary snarls file -> json snarls file
    with open(json_file, "w") as out:
        subprocess.run(
            ["vg", "view", "-R", str(snarls_file)],
            stdout=out,
            check=True,
            text=True,
        )


def compute_bubblegun_bubbles(tmpdir: Path):
    """
    Calls BubbleGun to detect bubbles and bubble chains
    and outputs the results in JSON format.
    """

    gfa_file = tmpdir / "compressed_graph.gfa"
    json_file = tmpdir / "bubbles.json"
    log_file = tmpdir / "bubblegun.log"

    with open(log_file, "w") as log:
        subprocess.run(
            [
                "BubbleGun",
                "-g", str(gfa_file),
                "bchains",
                "--bubble_json", str(json_file),
                "--save_memory",
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            text=True,
        )

    return json_file


def compute_panbubbles(tmpdir: Path):
    """
    Calls Billi to detect panbubbles and hairpins
    and outputs the results in a tab-separated text file.

    Parameters
    ----------
    tmpdir : Path
        Path to the directory containing the compressed graph.

    Returns
    -------
    output_file : Path
        Path to the file containing the detected panbubbles and hairpins.
    """

    gfa_file = (tmpdir / "compressed_graph.gfa").resolve()
    output_file = (tmpdir / "bubbles.json").resolve()
    log_file = (tmpdir / "billi.log").resolve()

    # Detect panbubbles and hairpins
    with open(output_file, "w") as out, open(log_file, "w") as log:
        subprocess.run(
            [
                "billi",
                "decompose",
                "-i", str(gfa_file),
                "-w",
            ],
            stdout=out,
            stderr=log,
            cwd=tmpdir.resolve(),
            check=True,
            text=True,
        )

    return output_file

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
                key = (str(end), str(start))
            else :
                key = (str(start), str(end))

            parent = None
            if "parent" in data:
                p = data["parent"]
                p_start = int(p["start"]["node_id"])
                p_end = int(p["end"]["node_id"])
                if p_start > p_end:
                    parent = (str(p_end), str(p_start))
                else :
                    parent = (str(p_start), str(p_end))

            boundaries.append(key)
            parents.append(parent)

    return boundaries, parents

def parse_bubblegun_bubbles(json_path):
    """
    Parse a BubbleGun JSON output file.

    Returns
    -------
    boundaries : list[tuple[str, str]]
        Sorted (source, sink) endpoint pairs for each bubble.

    parents : list[tuple[str, str] | None]
        Sorted (source, sink) endpoint pairs of the parent superbubble.
        None if the bubble has no parent.

    Both lists are aligned: boundaries[i] and parents[i]
    correspond to the same bubble.
    """

    with open(json_path, "r") as f:
        data = json.load(f)

    # Build a mapping from bubble IDs to their sorted endpoint pairs
    bubble_ends = {}

    for chain in data.values():
        for bubble in chain["bubbles"]:
            bubble_id = bubble["id"]
            start, end = bubble["ends"]
            if int(start) > int(end) :
                bubble_ends[bubble_id] = (end, start)
            else :
                bubble_ends[bubble_id] = (start, end)

    boundaries = []
    parents = []

    # Extract bubble boundaries and their parent superbubble boundaries
    for chain in data.values():
        parent_id = chain.get("parent_sb")

        if parent_id is not None:
            parent_ends = bubble_ends[parent_id]
        else:
            parent_ends = None

        for bubble in chain["bubbles"]:
            start, end = bubble["ends"]
            if int(start) > int(end) :
                boundaries.append((end, start))
            else :
                boundaries.append((start, end))
            parents.append(parent_ends)

    return boundaries, parents

def parse_panbubbles(file_path):
    """
    Parse a Billi output file.

    Returns
    -------
    boundaries : list[tuple[str, str]]
        Sorted endpoint pairs for each panbubble and hairpin.
        Hairpins are represented by (start, start).

    parents : list[tuple[str, str] | None]
        Sorted endpoint pairs of the parent bubble or hairpin.
        None if the bubble has no parent.

    Node IDs are returned as strings without orientation characters.
    Both lists are aligned.
    """

    def sorted_ends(ends):
        """Remove orientation characters and sort node IDs numerically."""
        node_ids = [
            node.lstrip("><")
            for node in ends
        ]
        return tuple(sorted(node_ids, key=int))

    records = []

    with open(file_path, "r") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")

            if not fields or fields[0] not in ("BB", "HP"):
                continue

            record_type = fields[0]

            if record_type == "BB":
                _, bubble_id, depth, parent_id, side1, side2, *_ = fields
                ends = (side1, side2)

            else:
                _, bubble_id, side1, side2, *_ = fields
                parent_id = None
                ends = (side1, side1)

            records.append({
                "type": record_type,
                "id": bubble_id,
                "parent_id": parent_id,
                "ends": sorted_ends(ends),
            })

    # Build a mapping from prefixed bubble IDs to their endpoint pairs
    bubble_ends = {
        f"{record['type']}:{record['id']}": record["ends"]
        for record in records
    }

    boundaries = []
    parents = []

    for record in records:
        boundaries.append(record["ends"])

        parent_id = record["parent_id"]

        if parent_id is None or parent_id == "-1":
            parents.append(None)
        else:
            parents.append(bubble_ends[parent_id])

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
    bubble_detector,
    save_subgraphs,
    compression,
    bubble_detector_timestamp,
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
    bubble_detector_time = bubble_detector_timestamp - compression["temporary_saving_time"]
    bubbles_load_and_index_time = bubbles["bubbles_load_and_index_time"] - bubble_detector_timestamp
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
        f.write(f"    Bubble detection tool             : {bubble_detector}\n")
        f.write(f"    Bubbles detected                  : {bubbles['nb_bubbles_detected']}\n")
        f.write(f"    Bubbles simplified                : {bubbles['nb_bubbles_simplified']}\n")
        f.write(f"    Removed nodes                     : {bubbles['removed_small_variants']}\n")
        f.write(f"    Remaining nodes                   : {bubbles['after_small_variants']}\n\n")

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
        f.write(f"Bubble detection time ({bubble_detector})     : {bubble_detector_time:.2f} s\n")
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
