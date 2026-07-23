"""
Here, we aim to extend the current GFA tag format by adding tags
that do respect the GFA naming convention.
A JSON string, PO (Path Offset) positions, relative to paths.
Hence, PO:J:{'w1':[(334,335,'+')],'w2':[(245,247,'-')]} tells that the walk/path w1
contains the sequence starting at position 334 and ending at position 335,
and the walk/path w2 contains the sequence starting at the offset 245 (ending 247),
and that the sequences are reversed one to each other.
Note that any non-referenced walk in this field means that the node
is not inside the given walk.
"""
import re
from gfagraphs import Graph
from collections import defaultdict
from pathlib import Path

from pgsimplify.utils import load_subgraphs

def compute_offsets(graph, subgraphs):
    """
    Adds to each node a JSON string containing the position of the begin and the end of the node on each path
    For example, the line S	872	CA PO:J:{'p1':[(238,240,'+')],'p2':[(312,314,'-')]} tells :
         - like in all graphs in gfa format, that node 872 contains sequence CA
         but also :
         - path p1 contains this sequence from position 238 to 240 
         - path p2 contains the reverse complement of this sequence from position 312 to 314
         - no other paths contains this node
    This function specially handle the case of simplified graph, using in addition to the graph a set of subgraphs, 
    corresponding to special nodes in the graph with name sgX with X an integer, for which the corresponding subgraph need to be traversed to compute real positions

    Parameters
    ----------
    graph : Graph
            main graph, containing abstracted node sgX with X an integer
    subgraphs : dict[str, Graph]
        subgraphs named sgX with X an integer corresponding to nodes of same name in the main graph
    """
    # Erase old positions in case they were some
    for seg_data in graph.segments.values():
        seg_data.pop("PO", None)

    # Pattern to recognize special nodes sgX
    sg_pattern = re.compile(r"^sg\d+$")

    # Compute positions for each path
    for path_name, path_data in graph.paths.items():
        current_pos = 0

        # Handle the case of a path traversing a subgraph multiple times
        sg_occurrences = defaultdict(int)

        # Compute position for each node of the path
        for node, orient in path_data["path"]:

            # Cases where the node is abstracted
            if sg_pattern.match(str(node)):
                
                # Check subgraph existence
                if node not in subgraphs:
                    raise ValueError(
                        f"Subgraph '{node}' referenced in path "
                        f"'{path_name}' but not found."
                    )

                # Compute path name
                sg = subgraphs[node]
                occ = sg_occurrences[(path_name, node)]
                sg_occurrences[(path_name, node)] += 1
                sg_path_name = f"{path_name}@{occ}"

                # Check path existence in the subgraph
                if sg_path_name not in sg.paths:
                    raise ValueError(
                        f"Path '{sg_path_name}' not found in subgraph '{node}'."
                    )

                # Compute size of the path in the graphs, equal size of the node
                node_length = sum(
                    sg.segments[sg_node]["length"]
                    for sg_node, _ in sg.paths[sg_path_name]["path"]
                    if sg_node not in {"source", "sink"}
                )

            # Case where the node is not abstracted   
            else:
                
                # Checks node existence in the graph
                if node not in graph.segments:
                    raise ValueError(
                        f"Node '{node}' referenced in path "
                        f"'{path_name}' but not found in graph."
                    )

                # Compute node size
                node_length = graph.segments[node]["length"]

            # Positions of the node on this path
            start = current_pos
            end = current_pos + node_length

            # Adds positions
            if "PO" not in graph.segments[node]:
                graph.segments[node]["PO"] = {}

            if path_name not in graph.segments[node]["PO"]:
                graph.segments[node]["PO"][path_name] = []

            graph.segments[node]["PO"][path_name].append(
                (
                    max(0, start),
                    end,
                    orient.value,
                )
            )

            current_pos = end

    # Ensure "PO" will not be recomputed by gfagraphs and will be saved
    graph.metadata["PO"] = True

def pipeline_offsets(input_dir : str, output_file: str):
    """
    Load data for computing offsets and saves the obtained graph

    Parameters
    ----------
    input_dir : Path
        Directory containing the files needed : main graph and subgraph directory 
    output_file : Path
        File to save the obtained graph
    """
    input_dir = Path(input_dir)
    if input_dir.is_dir():
        print("Loading main graph...")
        main_graph_path = input_dir / "main_graph.gfa"
        if not main_graph_path.exists():
            raise FileNotFoundError(f"Missing required input: {main_graph_path}")
        graph = Graph(str(main_graph_path), with_sequence=True)

        print("Loading subgraphs...")
        subgraphs_dir = input_dir / "subgraphs"
        if not subgraphs_dir.exists():
            raise FileNotFoundError(f"Missing required input: {subgraphs_dir}")
        subgraphs = load_subgraphs(subgraphs_dir)
        print(f"Loaded {len(subgraphs)} subgraphs")

    
    elif input_dir.suffix == ".gfa":
        if not input_dir.exists():
            raise FileNotFoundError(f"Missing graph: {input_dir}")
        graph = Graph(str(input_dir), with_sequence=True)
        subgraphs = {}

    else :
        raise ValueError(
        "Input must be either an existing directory or a '.gfa' file."
    )

    print("Computing offsets...")
    compute_offsets(
        graph=graph,
        subgraphs=subgraphs,
    )

    print("Saving graph...")
    graph.save_graph(output_file)

    print("Done.")