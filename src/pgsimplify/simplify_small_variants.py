import os
from gfagraphs import Graph
from collections import defaultdict
from itertools import pairwise
from pgsimplify.snarl_tree import SnarlTree

def extract_subgraphs(graph, pairs):
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
    pairs : (str, str)
        couple of nodes that are a bubble's boundaries

    Returns
    -------
    dict(list(str))
        Dictionnary containing for each bubble its interior nodes
        The keys are the index of bubbles in the snarl_pairs list
    """

    # Create dictionnaries start/end boundaries -> index in snarl_pairs list
    start_map = defaultdict(list)
    end_map   = defaultdict(list)
    for i, (s, t) in enumerate(pairs):
        start_map[s].append(i)
        end_map[t].append(i)

    # Initialize strutures to save interior nodes and bubbles sizes
    results = [set() for _ in pairs]
    snarls_size = [0 for _ in pairs]
    
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
            } for pid in range(len(pairs))
        }
    
        # Course of the path
        for node, _ in path_data['path']:

            # Current node is added to all opened intervals and its length to current length count
            for pid in active:
                if state[pid]["mode"] != 0:
                    results[pid].add(node)
                    state[pid]["length"] += graph.segments[node]["length"]

            # If the node is a bubble source
            if node in start_map:
                # For all pairs starting by this node
                for pid in start_map[node]:
                    # Interval is opened if it was closed
                    if pid not in active:
                        active.add(pid)
                        state[pid]["mode"] = 1 
                        results[pid].add(node) 
                        state[pid]["length"] += graph.segments[node]["length"]
                        
                    # Interval is closed if it was open by sink
                    elif pid in active and state[pid]["mode"] == -1:
                        active.remove(pid)
                        snarls_size[pid] = max(snarls_size[pid], state[pid]["length"]) # Update bubble size for the maximum
                        state[pid]["length"] = 0
                        
                    # If it was open by source, interval is closed and size is set to infinity 
                    # so that it won't be simplified because it does not correspond to a variant definition
                    elif pid in active and state[pid]["mode"] == 1:
                        active.remove(pid)
                        snarls_size[pid] = float('inf')
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
                        results[pid].add(node)
                        state[pid]["length"] += graph.segments[node]["length"]

                     # Interval is closed if it was open by source
                    elif pid in active and state[pid]["mode"] == 1:
                        active.remove(pid)
                        snarls_size[pid] = max(snarls_size[pid], state[pid]["length"]) # Update bubble size for the maximum
                        state[pid]["length"] = 0

                    # If it was open by sink, interval is closed and size is set to infinity 
                    # so that it won't be simplified because it does not correspond to a variant definition
                    elif pid in active and state[pid]["mode"] == -1:
                        active.remove(pid)
                        snarls_size[pid] = float('inf')
                        state[pid]["length"] = 0
                        #print(f"Warning : Pendant le parcours du path {path_name}, l'intervalle du sous graphe {pid}, de borne sortante {node} a été ouvert deux fois sans être fermé")
                    else:
                        print(f"Warning : Unkown interval state : {state[pid]["mode"]} (should be 0, 1, -1)")

        # At the end of the path, if interval are still open, it means the bubble doesn't correspond to variant definition
        # Size is set to infinity so that the bubble will nto be simplified
        for pid in active:
            snarls_size[pid] = float('inf')
            state[pid]["length"] = 0
            #print(f"Warning : Pendant le parcours du path {path_name}, l'intervalle du sous graphe {pid}, de borne sortante {node} n'a jamais été fermé")


    # Remove boudaries size from the bubble size
    for pid, (source, sink) in enumerate(pairs):
        snarls_size[pid] -= graph.segments[source]["length"]
        snarls_size[pid] -= graph.segments[sink]["length"]
        if snarls_size[pid] < 0:
            print("Warning : Negative bubble size")

    return results, snarls_size


def export_subgraphs(graph: Graph, chains, out_dir: str) -> None:
    """
    Export subgraphs using given bubble chains

    Parameters
    ----------
    graph : Graph
        Original graph
    chains : list[list[str, set(str), str]]
        Bubble chains, where each element has the form (source, interior nodes, sink)
    out_dir : str
        Output directory
    """
    # Create mapping : node -> id of the subgraph it is part of 
    subgraph_nodes = []
    node_to_subgraphs = defaultdict(list)

    for i, chain in enumerate(chains):
        nodes = set()

        for source, interior, sink in chain:
            nodes.add(source)
            nodes.add(sink)
            nodes.update(interior)

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

        nodes_in_path = {n for n, _ in pdata["path"]}

        candidate_sgs = set()

        for n in nodes_in_path:
            candidate_sgs.update(node_to_subgraphs.get(n, []))

        for sid in candidate_sgs:

            node_set = subgraph_nodes[sid]

            chunks = []
            current_chunk = []

            for n, o in pdata["path"]:

                if n in node_set:
                    current_chunk.append((n, o))

                else:
                    if current_chunk:
                        chunks.append(current_chunk)
                        current_chunk = []

            if current_chunk:
                chunks.append(current_chunk)

            # Introduces occurences in cases a same path traverse a subgrpahs multiple times
            for occ, chunk in enumerate(chunks):

                if not chunk:
                    continue

                subgraphs[sid].paths[f"{path_name}@{occ}"] = {
                    **pdata,
                    "path": chunk
                }

    # Save subgraphs
    for sid, sg in enumerate(subgraphs):

        sg.metadata = graph.metadata
        sg.headers = graph.headers
        output_path = os.path.join(out_dir, f"sg{sid}.gfa")
        sg.save_graph(output_path, minimal=True)


def build_bubbles_chains(bubbles):
    """
    Build chains of bubbles when boudaries overlap

    Parameters
    ----------
    bubbles : list[(str, set(str), str)]
        bubbles that will become element of a chain
    Returns
    -------
    list[list[(str, set(str), str)]]    
        the computed bubbles chains
    """
    # Build map : node -> number of time being a sink
    sink_count = defaultdict(int)
    for source, _ , sink in bubbles:
        sink_count[sink] += 1

    # Build map : source -> associated bubble
    source_to_bubble = {}
    starts = []

    for bubble in bubbles:
        source, _, sink = bubble
        source_to_bubble[source] = bubble
        if sink_count[source] == 0:
            starts.append(bubble) 

    # Build chains
    chains = []
    for start in starts:
        chain = [start]
        current = start
        while True:
            source, _, sink = current
            nxt = source_to_bubble.get(sink)
            if nxt is None:
                break
            chain.append(nxt)
            current = nxt
        chains.append(chain)

    return chains
        

def abstract_snarls(graph, snarl_nodes: list[set[str]], snarl_pairs, snarl_sizes, output_dir, save_subgraphs=True):
    """
    Compute simplified graph 

    Parameters
    ----------
    graph : Graph
        Graph to simplify
    snarl_nodes : list[set[str]]
        Nodes that are part of bubbles
    snarl_pairs : list[(str,str)]
        Nodes that are bubbles boundaries
    output_dir : str
        Output directory to store simplified graph and subgraphs
    save_subgraphs : bool
        Option to choose to save the subgraphs or not
    """
    # Changes the bubbles structure
    snarls_to_id = {}
    snarls = []
    i = 0
    for (source, sink), nodes in zip(snarl_pairs, snarl_nodes):
        internal_nodes = frozenset(nodes - {source, sink})
        snarl = (source, internal_nodes, sink)
        snarls.append(snarl)
        snarls_to_id[snarl] = i
        i += 1

    # Build chains
    chains = build_bubbles_chains(snarls)

    # Export of the subgraphs if the option is activated
    print("Save subgraphs...")
    if save_subgraphs :
        output_sg = os.path.join(output_dir, f"subgraphs")
        export_subgraphs(graph, chains, output_sg)

    # Compute simplified graph
    print("Compute simplified graph...")
    suppressed_nodes = set().union(*snarl_nodes)
    snarls_map = defaultdict(list)
    for i in range(len(chains)):
        chain = chains[i]
        # Create new node and map it to old one
        sg = f"sg{i}"
        snarls_map[chain[0][0]] = sg
        new_seq = graph.segments[chain[0][0]]["seq"]
        for snarl in chain:
            new_seq = new_seq + "N" * snarl_sizes[snarls_to_id[snarl]] + graph.segments[snarl[2]]["seq"]
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

        for n, o in pdata['path']:
            if n in snarls_map :
                new_path.append((snarls_map[n],o))
                continue
            if n in suppressed_nodes:
                continue

            new_path.append((n, o))

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
    

def compute_edge_orientation(graph):

    edge_orient = defaultdict(set)

    for path in graph.paths.values():

        for (u, ou), (v, ov) in pairwise(path["path"]):

            edge_orient[(u, v)].add(
                (ou.value, ov.value)
            )

    return edge_orient

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

    # Build snarl tree
    snarl_tree = SnarlTree(snarl_tree_path)
    snarl_pairs = list(snarl_tree.index_to_key.values())
    snarl_pairs = [(str(s), str(t)) for s, t in snarl_pairs]

    # Case where there is no snarl to simplify
    if snarl_pairs == []:
        print("No variants to simplify")
        return

    # Load graph
    print("Load graph...")
    graph = Graph(gfa_path)

    total_node_number = len(graph.segments)

    # Compute interior nodes and size of bubbles
    snarl_nodes, snarl_sizes = extract_subgraphs(graph, snarl_pairs)

    # Supression mask associated to bubbles index
    suppression_mask = [True if size > 0 and size < min_variant_size else False for size in snarl_sizes]
    print("Select variants to simplify...")
    # DFS of bubble tree
    snarls_to_abstract_index=[]
    stack = [snarl_tree.root]
    while stack:

        node = stack.pop()
        if node.key == (-1, -1): # artificial root of the tree
            for child in reversed(node.children):
                stack.append(child)
            continue

        i = snarl_tree.key_to_index[node.key]
        if suppression_mask[i]:
            snarls_to_abstract_index.append(i)
            skip =  True # subbubbles are skipped if the bubbles is already removed from the graph
        else :
            skip =  False

        if skip:
            continue
        # Reverse to keep natural order
        for child in reversed(node.children):
            stack.append(child)

    snarls_to_abstract = [nodes for i, nodes in enumerate(snarl_nodes) if i in snarls_to_abstract_index]
    snarl_pairs_to_abstract = [pairs for i, pairs in enumerate(snarl_pairs) if i in snarls_to_abstract_index]
    snarl_sizes = [sizes for i, sizes in enumerate(snarl_sizes) if i in snarls_to_abstract_index]

    # Abstract snarls in original graph and store subgraphs
    abstract_snarls(graph, snarls_to_abstract, snarl_pairs_to_abstract, snarl_sizes, output_dir, save_subgraphs)

    # Save final graph
    print(f"Save simplified graph...")
    os.makedirs(output_dir, exist_ok=True)
    output_gfa = os.path.join(output_dir, f"main_graph.gfa")
    graph.save_graph(output_gfa, minimal=True)

    # Simplification summary
    total_node_number_end = len(graph.segments)
    print(f"Small variant simplification : removed {total_node_number - total_node_number_end} nodes from the graph.")
    print(f"Number of nodes after simplification : {total_node_number_end}")
    return total_node_number_end

