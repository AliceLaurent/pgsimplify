from gfagraphs import Graph
from itertools import pairwise
from collections import defaultdict
from pgsimplify.utils import compute_edge_orientation, revcomp, get_substitutor


def is_same_type(type, succ, source, orient):
    """
    Verify if the edge having orientation orient is compatible with the existence of the branch of SNP/MNP described by type, succ and source
    An edge is not compatible, or the same type of an SNP/MNP branch if its existence create a conflict in the compression of the SNP/MNP

    Parameter
    ----------
    orient : 
        Orientation of the edge to compare
    type : type of the SNP/MNP branch, True for direct and False for indirect
    succ : True if the branch is a successor and false if it's a predecessor
    source : True if the branch is on the source side and sink if it's on the sink side

    Returns
    -------
    bool
        returns True if both edges are the same time, else returns False

    """
    if len(orient) == 0:
        return True
    if (succ and source and type == 'direct') or \
       (source and not succ and type == 'indirect') or \
       (not source and not succ and type == 'direct') or \
       (not source and succ and type ==  'indirect'):
        if ('+','+') in orient :
            return True
    else :
        if ('-','-') in orient :
            return True
    if (type == 'direct' and source) or (not source and type == 'indirect') :
        if ('+','-') in orient:
            return True
    else:
        if ('-','+') in orient :
            return True
    return False



def get_bubbles(graph: Graph, max_len=50):
    """
    Find all SNP/MNP structures in the graph (MNP under size max_len)

    Parameter
    ----------
    graph : Graph
        The pangenome graph to simplify
    max_len : int
        The maximum length of the MNP branches, by default 1 (means SNP detection only)

    Returns
    -------
    List[(str, List[str], str)]
        List of all SNP/MNP structures of the graph
    """

    # Compute for each edge both orientations of nodes at extremity of it
    edge_orientations = compute_edge_orientation(graph)

    # Init bubble list
    bubbles = []

    # Prevent multiple detections
    seen = set()

    # Nodes course
    for source in graph.segments:
        
        if source in seen:
            continue

        branches_type = {} # direct in cases of positive successors or negative predecessors, else indirect
        branches_orient = {} # Orientation of nodes at extremity of branches
        reversed_branches = {} # in cases source and sink need to be switched
        is_bubble_possible = True # becomes false if a condition for a bubble is not respected

        # Compute and course of neighbours that are successors
        neighbour_succ = set(graph.segments[source].get("successors", []))
        for n in neighbour_succ:
            branches_orient[n] = edge_orientations[(source,n)]
            if ('+','+') in branches_orient[n] or ('+','-') in branches_orient[n]:
                branches_type[n] = 'direct'
            else :
                branches_type[n] = 'indirect'

            if ('-','+') in branches_orient[n] or ('+','-') in branches_orient[n]:
                reversed_branches[n] = True
            else :
                reversed_branches[n] = False

         # Compute and course of neighbours that are predecessors
        neighbour_pred = set(graph.segments[source].get("predecessors", []))
        for n in neighbour_pred:
            branches_orient[n] = edge_orientations[(n,source)]
            if ('-','-') in branches_orient[n] or ('+','-') in branches_orient[n]:
                branches_type[n] = 'direct'
            else :
                branches_type[n] = 'indirect'
            
            if ('-','+') in branches_orient[n] or ('+','-') in branches_orient[n]:
                reversed_branches[n] = True
            else :
                reversed_branches[n] = False

        # Neighbour of branches are potential sinks
        neighbours_of_neighbours = defaultdict(set)

        # Compute neighbour of branches that are successors
        for branch in neighbour_succ:
            neighbours_of_neighbours[branch] = set(graph.segments[branch].get("successors", []))

         # Compute neighbour of branches that are predecessors
        for branch in neighbour_pred:
            neighbours_of_neighbours[branch] = set(graph.segments[branch].get("predecessors", [])) 
        
        # Create dictionnary : branches -> potential sink 
        potential_sinks = defaultdict(set)
        for branch, neighbours in neighbours_of_neighbours.items():
            for neighbour_of_neighbour in neighbours:
                potential_sinks[neighbour_of_neighbour].add(branch)
        
        if not is_bubble_possible:
            continue

        # Course of all potential sinks
        for sink, branches in potential_sinks.items():
            
            if sink in seen: # May have been the sink of another SNP/MNP
                continue

            is_bubble_possible = True # each sink is independant

            if source == sink: # source and sink have to be different
                continue

            if len(branches) < 2: # a SNP/MNP has at least to branches
                continue

            # Verify each branch has at least 2 neighbours
            # Compute branches lengths
            # Verify that all branches have the same orientation (needed to form a coherent SNP/MNP)
            lengths = set()
            bubble_orient = None
            for branch in branches:
                neigh = set(graph.segments[branch].get("successors", [])) | set(graph.segments[branch].get("predecessors", []))
                if len(neigh) > 2:
                    is_bubble_possible = False
                    break
                if bubble_orient is None :
                    bubble_orient = branches_type[branch]
                else : 
                    if branches_type[branch] != bubble_orient:
                        is_bubble_possible = False
                        break

                branch_len = graph.segments[branch]["length"]
                lengths.add(branch_len)
            
            if is_bubble_possible : 
                # Verify all branches have same length
                if len(lengths) > 1:
                    continue

                # Verify branches length is inferior to defined max_len
                if next(iter(lengths)) > max_len:
                    continue

                # Verify that no successor edge of source makes SNP/MNP compression unpossible
                source_succ = set(graph.segments[source].get("successors", []))
                for succ in source_succ:
                    if succ in branches:
                        continue
                    succ_orient = edge_orientations[(source,succ)]
                    if is_same_type(branches_type[next(iter(branches))], True, True, succ_orient):
                        is_bubble_possible = False

                # Verify that no predecessor edge of source makes SNP/MNP compression unpossible
                source_pred = set(graph.segments[source].get("predecessors", []))
                for pred in source_pred:
                    if pred in branches:
                        continue
                    pred_orient = edge_orientations[(pred,source)]
                    if is_same_type(branches_type[next(iter(branches))], False, True, pred_orient):
                        is_bubble_possible = False

                # Verify that no successor edge of sink makes SNP/MNP compression unpossible
                sink_succ = set(graph.segments[sink].get("successors", []))
                for succ in sink_succ:
                    if succ in branches:
                        continue
                    succ_orient = edge_orientations[(sink,succ)]
                    if is_same_type(branches_type[next(iter(branches))], True, False, succ_orient):
                        is_bubble_possible = False

                # Verify that no predecessor edge of sink makes SNP/MNP compression unpossible
                sink_pred = set(graph.segments[sink].get("predecessors", []))
                for pred in sink_pred:
                    if pred in branches:
                        continue
                    pred_orient = edge_orientations[(pred,sink)]
                    if is_same_type(branches_type[next(iter(branches))], False, False, pred_orient):
                        is_bubble_possible = False
                
                if is_bubble_possible:
                    # Compute revcomp of branch sequence if it' reversed
                    for branch in branches :
                        if reversed_branches[branch]:
                            graph.segments[branch]['seq'] = revcomp(graph.segments[branch]['seq'])
                    # Order source and sink depending on branches type direct or indirect
                    if branches_type[next(iter(branches))] == 'direct':
                        bubbles.append((source, branches, sink))
                    else :
                        bubbles.append((sink, branches, source))
                    
                    # Mark source and branches as seen
                    seen.add(source)
                    for branch in branches:
                        seen.add(branch)
            
    return bubbles


def build_bubbles_chains(bubbles):
    """
    Form SNP/MNP chains
    When a node is at the same time source of a SNP/MNP and sink of another, both structures can be chained together

    Parameter
    ----------
    bubbles : 
        List of SNP/MNP structures with form (source : int, branches : int list, sink : int)
    Returns
    -------
    List[List[(str, List[str], str)]]
        List of bubbles organized in chains
    """
    # Create dictionnary for number of time a node appears as sink in a SNP/MNP structure
    sink_count = defaultdict(int)
    for source, _ , sink in bubbles:
        sink_count[sink] += 1

    # Create dictionnary source -> bubble and a list of all chains starts (nodes that are not sink in any structure)
    source_to_bubble = {}
    starts = []

    for bubble in bubbles:
        source, _, sink = bubble
        source_to_bubble[source] = bubble
        if sink_count[source] == 0:
            starts.append(bubble) 

    # Build chain usinge source -> bubble dictionnary and start list
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


def compress_bubbles_chains(graph, max_len=1):
    """
    Compress graph by replacing SNP/MNP chains by a single compressed node

    Parameter
    ----------
    graph : Graph
        The pangenome graph to compress
    max_len : the maximum length of a MNP to compress

    Returns
    -------
    dict(int)
        Offsets dictionnary
    """
   
    # Compute SNP/MNP structures
    bubbles = get_bubbles(graph, max_len)

    # Organize structures in SNP/MNP chains
    chains = build_bubbles_chains(bubbles)

    # For cases where a path begins in the middle of a future compressed node
    offsets = {}
    first_node_to_paths = defaultdict(list)
    first_nodes = {}

    # Build dictionnary first node of a path -> path name
    for name, pdata in graph.paths.items():
        if pdata['path']:
            first_node = pdata['path'][0][0]
            first_node_to_paths[first_node].append(name)

    # Compute nodes that will be removed and new sequence that will replace them for each chain
    suppressed_nodes = set()
    for chain in chains:

        # First node of the chain will replace the entire chain
        final_node,_,_ = chain[0]
        # It will receive a new sequence replacing all nodes of the chain
        new_seq = graph.segments[final_node]['seq']

        # List for the path having their first node in this chain
        offset_paths = []

        # Course of the chain
        for  source, branches, sink in chain:
            
            # Compute sequence to replace the chain
            seqs = []
            for b in branches:
                s = graph.segments[b]['seq']
                seqs.append(s)

            seqs = list(set(seqs))
            sub = get_substitutor(seqs)

            # Add chain element to new sequence
            new_seq = new_seq + sub + graph.segments[sink]['seq']

            # Add sink and branches to supressed nodes and verify if any node of this chain element is the beginning of a path
            suppressed_nodes.add(sink)

            for b in branches :
                suppressed_nodes.add(b)
                if b in first_node_to_paths :
                    offset_paths = first_node_to_paths[b]
                
            if source in first_node_to_paths :
                offset_paths = first_node_to_paths[source]
            if sink in first_node_to_paths :
                offset_paths = first_node_to_paths[sink]

        # In cases of an offset (path begin in the chain), a set of all nodes in the chain is computed
        chain_nodes = set()

        for source, branches, sink in chain:
            chain_nodes.add(source)
            chain_nodes.add(sink)

            for branch in branches:
                chain_nodes.add(branch)

        # Course of begin of paths that need an offset for this chain to count which part of the chain is run through by the path
        for path_name in offset_paths:
            len_of_chain_in_path = 0
            final_node_in_hap = False
            for n, o in graph.paths[path_name]['path']:
                if n not in chain_nodes:
                    break
                else :
                    len_of_chain_in_path += graph.segments[n]['length']
                    if n == final_node :
                        final_node_in_hap = True
            # Offset is substraction of length of the chain in path to length of compressed sequence
            len_offset = len(new_seq) - len_of_chain_in_path
            if len_offset != 0 :
                offsets[path_name] = len_offset
            # In cases of offset where the first node of the chain is not in the path
            # The first node of the chain is saved as first node of the path (since onmy the first node of the chain is saved in final graph)
            # Orientation is the one of the real first node of the path
            if not final_node_in_hap:
                first_nodes[path_name] = (final_node, graph.paths[path_name]['path'][0][1])
        
        # Content of the final is replaced by the new compressed content
        graph.segments[final_node]['seq'] = new_seq
        graph.segments[final_node]['length'] = len(new_seq)

    # Compute new segment set
    new_segments = {}

    for n, data in graph.segments.items():

        if n in suppressed_nodes : # Nodes of compressed chain (except first node of chain) are suppressed
            continue
        
        new_segments[n] = data # Other nodes are kept

    graph.segments = new_segments

    # Paths rewriting
    new_paths = {}

    for name, pdata in graph.paths.items():

        new_path = []
        # First node is artificially added if needed
        if name in first_nodes:
            new_path.append(first_nodes[name])

        # Path are copied except for the supressed nodes
        for n, o in pdata['path']:

            if n not in graph.segments:
                continue

            new_path.append((n, o))

        new_paths[name] = {
        **pdata,
        'path': new_path
    }
    graph.paths = new_paths

    # New edges set
    graph.lines = {}
    # Edges are built using paths
    for pdata in graph.paths.values():
        for (u, ou), (v, ov) in pairwise(pdata["path"]):
            graph.add_edge(u, ou, v, ov)


    # Update of sequence offsets and nieghbours
    graph.sequence_offsets()
    graph.compute_neighbors()

    return offsets


def iterative_bubble_compression(graph, offsets, max_len_to_collapse = 50, nb_iter = 2):
    """
    Iterate SNP/MNP Compression
    
    Parameter
    ----------
    graph : Graph
        Graph to compress
    offsets : 
        Offsets dictionnary
    nb_iter : int
        Number of iteration for compression
    max_len_to_collapse : int
        Maximum length for MNP simplification

    Returns
    -------
    dict(int)
        Offsets dictionnary
    """
    total_node_number_before_iter = len(graph.segments)
    for i in range(nb_iter):
        print(f"--- Compression iteration {i+1} ---")
        
        # Call compression 
        offset_iter = compress_bubbles_chains(
            graph=graph,
            max_len=max_len_to_collapse
        )

        # Add new offsets
        new_offsets = offsets.copy()

        for path, value in offset_iter.items():
            new_offsets[path] = new_offsets.get(path, 0) + value

        offsets = new_offsets

        # Print iteration compression summary
        nb_segments = len(graph.segments)
        print(f"SNP/MNP compression : removed {total_node_number_before_iter - nb_segments} nodes form the graph")
        total_node_number_before_iter = nb_segments
    
    return offsets

