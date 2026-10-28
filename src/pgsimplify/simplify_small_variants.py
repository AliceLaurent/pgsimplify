import os
import time
from gfagraphs import Graph
from collections import defaultdict
from itertools import pairwise
from pgsimplify.snarl_tree import SnarlTree
from pgsimplify.utils import parse_snarls
from pgsimplify.bubble import Bubble

def extract_subgraphs(graph, bubbles, boundary_to_bubble_index):
    """
    Compute size and extact the nodes that are inside of a bubble represented by a couple of boundary nodes 
    For every buble, nodes that are on a path between the two boundary nodes will be added to its internal node sets
    To do so, path are traversed and each bubble has its own interval that is opened when the first boudary node is met and
    closed when the second is found. While the interval remains open, current nodes are added to internal nodes set
    The longest path between the two boudaries in term of number of bases found in the graph is the bubble size

    Parameters
    ----------
    graph : Graph
        Graph containing the bubbles
    bubbles : list(Bubble)
        List of bubbles
    boundary_to_bubble_index : dict((str, str) -> int)
        Mapping from boundary nodes to bubble index

    Returns
    -------
    dict(list(str))
        Dictionnary containing for each bubble its interior nodes
        The keys are the index of bubbles in the snarl_pairs list
    """
    # Create dictionnaries start/end boundaries -> index in snarl_pairs list
    start_map = defaultdict(list)
    end_map   = defaultdict(list)
    for i in range(len(bubbles)):
        s, t = bubbles[i].boundaries
        start_map[s].append(i)
        end_map[t].append(i)

    # Traverse all paths only once
    for _, path_data in graph.paths.items():
         # consigne l'etat d'un intervalle
        # 0 -> L'intervalle est fermé (aucune extremit n'a ete rencontre dans le path)
        # 1 -> On a  rencontre source mais pas sink (on est entree dans l'intervalle dans le "bon" sens)
        # -1 -> On a rencontre sink, mais pas source (on est entree dans le sous graphe par sink)
        # Active is the set of interval that are currently open
        active = set()

        # State is a dictionnary that store for each bubble :
        #   - mode : the direction the bubble was entered (1 for start met first, -1 for end met first, 0 if outside of the bubble)
        #   - length : used to count the size of a bubble on a path
        state = {
            pid: {
                "mode": 0,
                "length": 0
            } for pid in range(len(bubbles))
        }
    
        # Course of the path
        for node, orient in path_data['path']:

            # Current node is added to all opened intervals and its length to current length count
            for pid in active:
                if state[pid]["mode"] != 0:
                    bubbles[pid].interior_nodes.add(node)
                    state[pid]["length"] += graph.segments[node]["length"]

            # If the node is a bubble source
            if node in start_map:
                # For all pairs starting by this node
                for pid in start_map[node]:
                    # Interval is opened if it was closed
                    if pid not in active:
                        active.add(pid)
                        state[pid]["mode"] = 1 
                        bubbles[pid].interior_nodes.add(node)
                        state[pid]["length"] += graph.segments[node]["length"]
                        
                    # Interval is closed if it was open by sink
                    elif pid in active and state[pid]["mode"] == -1:
                        active.remove(pid)
                        bubbles[pid].size = max(bubbles[pid].size, state[pid]["length"]) # Update bubble size for the maximum
                        state[pid]["length"] = 0
                        
                    # If it was open by source, interval is closed and size is set to infinity 
                    # so that it won't be simplified because it does not correspond to a variant definition
                    elif pid in active and state[pid]["mode"] == 1:
                        active.remove(pid)
                        bubbles[pid].size = float('inf')
                        state[pid]["length"] = 0
                        #print(f"Warning : Pendant le parcours du path {path_name}, l'intervalle du sous graphe {pid}, de borne entrante {node} a été ouvert deux fois sans être fermé")
                    else:
                        print(f"Warning : Unkown interval state : {state[pid]["mode"]} (should be 0, 1, -1)")

            # If the node is a bubble sink
            if node in end_map:
                # For all pairs ending by this node
                for pid in end_map[node]:
                    # Interval is open if it was closed
                    if pid not in active:
                        active.add(pid)
                        state[pid]["mode"] = -1
                        bubbles[pid].interior_nodes.add(node)
                        state[pid]["length"] += graph.segments[node]["length"]

                     # Interval is closed if it was open by source
                    elif pid in active and state[pid]["mode"] == 1:
                        active.remove(pid)
                        bubbles[pid].size = max(bubbles[pid].size, state[pid]["length"]) # Update bubble size for the maximum
                        state[pid]["length"] = 0

                    # If it was open by sink, interval is closed and size is set to infinity 
                    # so that it won't be simplified because it does not correspond to a variant definition
                    elif pid in active and state[pid]["mode"] == -1:
                        active.remove(pid)
                        bubbles[pid].size = float('inf')
                        state[pid]["length"] = 0
                        #print(f"Warning : Pendant le parcours du path {path_name}, l'intervalle du sous graphe {pid}, de borne sortante {node} a été ouvert deux fois sans être fermé")
                    else:
                        print(f"Warning : Unkown interval state : {state[pid]["mode"]} (should be 0, 1, -1)")

        # At the end of the path, if interval are still open, it means the bubble doesn't correspond to variant definition
        # Size is set to infinity so that the bubble will not be simplified
        for pid in active:
            bubbles[pid].size = float('inf')
            state[pid]["length"] = 0
            #print(f"Warning : Pendant le parcours du path {path_name}, l'intervalle du sous graphe {pid}, de borne sortante {node} n'a jamais été fermé")


    # Remove boudaries size from the bubble size
    for pid in range(len(bubbles)):
        source, sink = bubbles[pid].boundaries
        bubbles[pid].size -= graph.segments[source]["length"]
        bubbles[pid].size -= graph.segments[sink]["length"]
        if bubbles[pid].size < 0:
            print("Warning : Negative bubble size")


def export_subgraphs(graph: Graph, chains, bubbles: dict[int, Bubble], out_dir: str) -> None:
    """
    Export subgraphs using given bubble chains

    Parameters
    ----------
    graph : Graph
        Original graph
    chains : list[list[int]]
        Bubble chains, where each element is a list of indices of bubbles in the original bubbles list
    bubbles : dict[int, Bubble]
        Dictionary of bubbles to simplify
    out_dir : str
        Output directory
    """
    # Create mapping : node -> id of the subgraph it is part of 
    subgraph_nodes = []
    node_to_subgraphs = defaultdict(list)

    for i, chain in enumerate(chains):
        nodes = set()

        for pid in chain:
            bubble = bubbles[pid]
            nodes.add(bubble.boundaries[0])
            nodes.add(bubble.boundaries[1])
            nodes.update(bubble.interior_nodes)

        subgraph_nodes.append(nodes)

        for n in nodes:
            node_to_subgraphs[n].append(i)

    # Initialize empty subgraphs
    subgraphs = []

    for _ in chains:
        sg = Graph()
        sg.segments = {}
        sg.lines = {}
        sg.paths = {}
        subgraphs.append(sg)

    # Add nodes to subgraphs
    for node, data in graph.segments.items():

        for sid in node_to_subgraphs.get(node, []):
            subgraphs[sid].segments[node] = data

    # Add edges to subgraphs
    for (u, v), data in graph.lines.items():

        common = (
            set(node_to_subgraphs.get(u, []))
            &
            set(node_to_subgraphs.get(v, []))
        )

        for sid in common:
            subgraphs[sid].lines[(u, v)] = data

    # Add paths to subgraphs
    for path_name, pdata in graph.paths.items():
        candidate_chunks = defaultdict(list)
        current_chunks = {}

        for n, o in pdata["path"]:
            subgraph_ids = node_to_subgraphs.get(n, [])
            active_ids = set(subgraph_ids)

            # Close chunks for subgraphs that do not contain this node
            for sid in list(current_chunks):
                if sid not in active_ids:
                    candidate_chunks[sid].append(current_chunks.pop(sid))

            # Extend or start a chunk for each subgraph containing this node
            for sid in subgraph_ids:
                if sid not in current_chunks:
                    current_chunks[sid] = []

                current_chunks[sid].append((n, o))

        # Close remaining chunks
        for sid, chunk in current_chunks.items():
            candidate_chunks[sid].append(chunk)

        # Store path occurrences, preserving the original naming convention
        for sid, chunks in candidate_chunks.items():
            for occ, chunk in enumerate(chunks):
                if not chunk:
                    continue

                subgraphs[sid].paths[f"{path_name}@{occ}"] = {
                    **pdata,
                    "path": chunk,
                }

    subgraph_build_time = time.perf_counter()

    # Save subgraphs
    for sid, sg in enumerate(subgraphs):

        sg.metadata = graph.metadata
        sg.headers = graph.headers
        output_path = os.path.join(out_dir, f"sg{sid}.gfa")
        sg.save_graph(output_path, minimal=True)
    subgraph_save_time = time.perf_counter()
    return subgraph_build_time, subgraph_save_time


def build_bubbles_chains(bubbles):
    """
    Build chains of bubbles when boundaries overlap.

    Parameters
    ----------
    bubbles : dict[int, Bubble]
        Bubbles with boundaries = (source, sink).

    Returns
    -------
    list[list[int]]
        Chains represented by indices in the bubbles list.
    """
    # Build map: endpoint -> index of the bubble
    endpoint_to_index = defaultdict(list)
    endpoint_count = defaultdict(int)

    for i, bubble in bubbles.items():
        source, sink = bubble.boundaries
        endpoint_to_index[source].append(i)
        endpoint_to_index[sink].append(i)
        endpoint_count[source] += 1
        endpoint_count[sink] += 1

    # Find chains endpoints (endpoints appearing only once in the bubbles)
    chain_endpoints = []
    for endpoint, count in endpoint_count.items():
        if count == 1:
            chain_endpoints.append(endpoint)
        if count > 2:
            raise ValueError(f"Node is boundary of more than two bubbles : {endpoint}")

    # Build chains
    chains = []
    visited_bubbles = set()

    for endpoint in chain_endpoints:
        start_bubble = endpoint_to_index[endpoint][0]
        if start_bubble in visited_bubbles:
            continue
        chain = [start_bubble]
        current_endpoint= endpoint
        current_bubble = start_bubble
        while True:
            visited_bubbles.add(current_bubble)
            source, sink = bubbles[current_bubble].boundaries
            if current_endpoint == source:
                neighbour_list = endpoint_to_index[sink]
            else:
                neighbour_list = endpoint_to_index[source]  
            if len(neighbour_list) < 2:
                break
            if neighbour_list[0] != current_bubble:
                nxt = neighbour_list[0]
            elif len(neighbour_list) > 1 and neighbour_list[1] != current_bubble:
                nxt = neighbour_list[1]
            else:
                raise ValueError(f"Cannot find next bubble in chain for bubble {current_bubble.boundaries} with sink {sink}")
            
            chain.append(nxt)
            current_bubble = nxt
            if current_endpoint == source:
                current_endpoint = sink
            else:
                current_endpoint = source

        chains.append(chain)

    missing = set(bubbles) - visited_bubbles
    if missing:
        raise ValueError(f"Bubbles not covered by any chain: {missing}")

    return chains
        

def abstract_snarls(graph, bubbles, output_dir, save_subgraphs=True):
    """
    Compute simplified graph 

    Parameters
    ----------
    graph : Graph
        Graph to simplify
    bubbles : dict[int, Bubble]
        Dictionary of bubbles to simplify, where the keys are the indices of the bubbles in the original bubbles list
    output_dir : str
        Output directory to store simplified graph and subgraphs
    save_subgraphs : bool
        Option to choose to save the subgraphs or not
    """

    # Build chains
    chains = build_bubbles_chains(bubbles)
    snarls_abstraction_time = time.perf_counter()

    # Export of the subgraphs if the option is activated
    subgraph_build_time, subgraph_save_time = 0, 0
    if save_subgraphs :
        print("Save subgraphs...")
        output_sg = os.path.join(output_dir, f"subgraphs")
        subgraph_build_time, subgraph_save_time = export_subgraphs(graph, chains, bubbles, output_sg)

    # Compute simplified graph
    print("Compute simplified graph...")
    node_to_sg = {}
    suppressed_nodes = set().union(*[bubbles[pid].interior_nodes for pid in bubbles])

    for i in range(len(chains)):
        chain = chains[i]
        # Create new node
        sg = f"sg{i}"

        # All snarl nodes are mapped to sg
        for i in range(len(chain)):
            pid = chain[i]
            bubble = bubbles[pid]
            for node in bubble.interior_nodes:
                node_to_sg[node] = sg

        new_seq = graph.segments[bubbles[chain[0]].boundaries[0]]["seq"]
        for pid in chain:
            bubble = bubbles[pid]
            new_seq = new_seq + "N" * bubble.size + graph.segments[bubble.boundaries[1]]["seq"]
        graph.segments[sg] = {
        "seq": new_seq,
        "length": len(new_seq),
        "successors": set(),
        "predecessors": set()
        }

    # Compute new node set
    new_segments = {}
    for n, data in graph.segments.items():
        if n in suppressed_nodes :
            continue
        new_segments[n] = data
    graph.segments = new_segments

    # Compute new paths
    new_paths = {}
    for name, pdata in graph.paths.items():

        new_path = []
        current_sg = None

        for n, o in pdata['path']:

            sg = node_to_sg.get(n)
            if sg is not None:
                # Replace only the first node met of the snarl
                if sg != current_sg:
                    new_path.append((sg, o))
                    current_sg = sg
            else:
                new_path.append((n, o))
                current_sg = None

        new_paths[name] = {
        **pdata,
        'path': new_path
    }
    graph.paths = new_paths

    # Compute new edges using paths
    graph.lines = {}

    for pdata in graph.paths.values():
        for (u, ou), (v, ov) in pairwise(pdata["path"]):
            graph.add_edge(u, ou, v, ov)

    # Update the graph
    graph.sequence_offsets()
    graph.compute_neighbors()
    main_graph_build_time = time.perf_counter()

    return snarls_abstraction_time, main_graph_build_time, subgraph_build_time, subgraph_save_time

def compress_snarls_pipeline(
    gfa_path: str,
    snarl_tree_path: str,
    output_dir: str,
    min_variant_size: int = 50,
    save_subgraphs: bool = True,
):
    """
    Extract small variants from the graph and store them as subgraphs
    1. Load bubbles and compute bubbles tree
    2. Load the graph
    3. Compute interior nodes and size of all bubbles
    4. Performs dfs in bubble tree to find the bubbles to simplify
    5. Organize the bubbles in chains, extract bubbles from original graph and save them
    6. Save the final graph

    Parameters
    ----------
    gfa_path : str
        input GFA path
    snarl_tree_path : str
        input bubble tree path
    output_gfa : str
        directory to store final simplified graph
    min_variant_size : int
        Variant smaller than min_variant_size will be simplified
    save_subgraphs : bool
        option tos ave or not the subgraphs
    """
    # Load bubbles
    boundaries_list, parents_list = parse_snarls(snarl_tree_path)

    # Initialize bubbles list and build hashtable to access bubbles indexes by their boundaries
    bubbles_list = []
    boundary_to_bubble_index = defaultdict(list)
    
    for i in range(len(boundaries_list)):
        bubbles_list.append(Bubble(boundaries_list[i]))
        start, end = bubbles_list[i].boundaries
        boundary_to_bubble_index[start].append(i)
        boundary_to_bubble_index[end].append(i)

    # Adds parent id to each bubble structure
    for i in range(len(bubbles_list)):
        parent_boundaries = parents_list[i]
        if parent_boundaries is None:
            continue
        start, end = parent_boundaries
        start_bubbles = boundary_to_bubble_index[start]
        end_bubbles = boundary_to_bubble_index[end]
        parent_indices = set(start_bubbles) & set(end_bubbles)

        # The parent of a bubble is unique, if there is no parent or multiple parents, an error is raised
        if len(parent_indices) != 1:
            raise ValueError(
                f"Parent not unique or not found for bubble {i}: "
                f"{parent_boundaries}"
            )

        parent_index = parent_indices.pop()
        bubbles_list[i].parent = parent_index


    # Case where there is no snarl to simplify
    if bubbles_list == []:
        print("No variants to simplify")
        return
    
    bubble_load_and_index_time = time.perf_counter()

    # Load graph
    print("Load graph...")
    graph = Graph(gfa_path)
    graph_load_time = time.perf_counter()

    total_node_number = len(graph.segments)

    # Compute interior nodes and size of bubbles
    extract_subgraphs(graph, bubbles_list, boundary_to_bubble_index)
    node_and_size_extraction_time = time.perf_counter()

    # Select the bubbles to simplify, which are the ones that are smaller than min_variant_size 
    # and whose parent not simplified (if the parent is simplified, it will contain the current bubble, so non need to simplify it)
    print("Select variants to simplify...")
    bubbles_to_abstract = {}
    for i, bubble in enumerate(bubbles_list):
        if bubble.size > 0 and bubble.size < min_variant_size :
            if bubble.parent is None :
                bubbles_to_abstract[i] = bubble
            else :
                parent = bubbles_list[bubble.parent] 
                if parent.size >= min_variant_size :
                    bubbles_to_abstract[i] = bubble
    
    bubbles_selection_time = time.perf_counter()

    # Abstract snarls in original graph and store subgraphs
    snarls_abstraction_time, main_graph_build_time, subgraph_build_time, subgraph_save_time = abstract_snarls(graph, bubbles_to_abstract, output_dir, save_subgraphs)

    # Save final graph
    print(f"Save simplified graph...")
    os.makedirs(output_dir, exist_ok=True)
    output_gfa = os.path.join(output_dir, f"main_graph.gfa")
    graph.save_graph(output_gfa, minimal=True)
    graph_save_time = time.perf_counter()

    # Simplification summary
    total_node_number_end = len(graph.segments)
    print(f"Small variant simplification : removed {total_node_number - total_node_number_end} nodes from the graph.")
    print(f"Number of nodes after simplification : {total_node_number_end}")
    
    return {
    "before_small_variants": total_node_number,
    "after_small_variants": total_node_number_end,
    "removed_small_variants": total_node_number - total_node_number_end,
    "nb_bubbles_detected": len(bubbles_list),
    "nb_bubbles_simplified": len(bubbles_to_abstract),
    "bubbles_load_and_index_time": bubble_load_and_index_time,
    "graph_load_time": graph_load_time,
    "node_and_size_extraction_time": node_and_size_extraction_time,
    "bubbles_selection_time": bubbles_selection_time,
    "bubbles_abstraction_time": snarls_abstraction_time,
    "main_graph_build_time": main_graph_build_time,
    "subgraph_build_time": subgraph_build_time,
    "subgraph_save_time": subgraph_save_time,
    "graph_save_time": graph_save_time
}

