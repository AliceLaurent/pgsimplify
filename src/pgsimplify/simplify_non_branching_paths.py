from gfagraphs import Graph
from itertools import pairwise
from collections import defaultdict

from pgsimplify.utils import compute_edge_orientation, invert_orient

def get_non_branching_nodes(graph : Graph):
    """
    Find nodes that could be part of a non-branching path

    Parameter
    ----------
    graph : Graph
        Graph to compress
    """
    edge_orient = compute_edge_orientation(graph)
    chain_nodes = []

    # Course of nodes
    for node in graph.segments:
        
        # Compute succs and preds of the node
        succs = set(graph.segments[node].get("successors", []))
        preds = set(graph.segments[node].get("predecessors", []))

        # Nodes can be begin, middle or end of a non-branching paths
        # It will be the middle if both sides are possible and have a neighbour candidate
        # It will be a begin if only right side is possible and an end if only left side is possible
        right_possible = True
        right_candidate = None
        left_possible = True
        left_candidate = None

        # Neighbours course
        for succ in succs:
            
            orient = edge_orient[(node,succ)]
            # Cases where only the equivalent edge is in the graph (opposite direction and opposite orientation)
            if len(orient) == 0:
                orient = edge_orient.get((succ, node), set())
                orient = invert_orient(orient)
            # Should not happen
            if len(orient) == 0:
                continue
            # May happen, but is not non-branching node
            if len(orient) != 1:
                right_possible = False
                left_possible = False
                #print(f"Warning : il y a plusieurs orientations entre {pred} et {node} : {orient}")
            if ('+','-') in orient:
                right_possible = False
            if ('-','+') in orient:
                left_possible = False
            if ('+','+') in orient:
                if right_candidate is not None :
                    right_possible = False
                else :
                    right_candidate = succ
            if ('-','-') in orient:
                if left_candidate is not None :
                    left_possible = False
                else :
                    left_candidate = succ
        
        if right_possible or left_possible:
            for pred in preds:
                orient = edge_orient[(pred,node)]
                # Cases where only the equivalent edge is in the graph (opposite direction and opposite orientation)
                if len(orient) == 0: 
                    orient = edge_orient.get((node, pred), set())
                    orient = invert_orient(orient)
                # Should not happen
                if len(orient) == 0:
                    continue
                # May happen, but is not non-branching node
                if len(orient) != 1:
                    right_possible = False
                    left_possible = False
                    #print(f"Warning : il y a plusieurs orientations entre {pred} et {node} : {orient}")
                if ('+','-') in orient:
                    right_possible = False
                if ('-','+') in orient:
                    left_possible = False
                if ('-','-') in orient:
                    if right_candidate is None :
                        right_candidate = pred
                    if right_candidate is not None and right_candidate != pred:
                        right_possible = False
                if ('+','+') in orient:
                    if left_candidate is None :
                        left_candidate = pred
                    elif left_candidate is not None and left_candidate != pred:
                        left_possible = False

        if left_possible and left_candidate is not None and right_possible and right_candidate is not None:
            chain_nodes.append((left_candidate, node, right_candidate))
        elif left_possible and left_candidate is not None :
            chain_nodes.append((left_candidate, node, None))
        elif right_possible and right_candidate is not None:
            chain_nodes.append ((None, node, right_candidate))

    return chain_nodes


def build_non_branching_paths(chain_nodes):
    """
    Build chains of non-branching nodes

    Parameter
    ----------
    chain_nodes : 
        (left, node, right) with node a non-branching node and left and right potential successors or predecessors in the chain
        Left and right can be None, and are linked to node but not necessarily meet the condition of non-branching node
       
    """
    # Sets of neighbour of non-branching nodes
    left_possible = set()
    right_possible = set()

    for left, node, right in chain_nodes:
        if left is not None :
            left_possible.add(node)
        if right is not None:
            right_possible.add(node)

    # Left and right nodes needs to be non-branching nodes in another structure, or it means that they don't meet the conditions
    # for a non-branching nodes and can't be part of a chain
    filtered_chain_nodes = []
    for left, node, right in chain_nodes:
        if (left is None or left not in right_possible) and right in left_possible:
            filtered_chain_nodes.append((None, node, right))
            continue
        
        if (right is None or right not in left_possible) and left in right_possible:
            filtered_chain_nodes.append((left, node, None))
            continue

        if left is not None and right is not None and right in left_possible and left in right_possible:
            filtered_chain_nodes.append((left,node,right))
            continue
        # Left and right both None means a chain of length 1, which can't be simplified

    # Dictionnary node -> left/right node
    left_of = {}
    right_of = {}

    for left, node, right in filtered_chain_nodes:

        left_of[node] = left
        right_of[node] = right

    visited = set()
    chains = []

    # Build chains
    for (left, node, right) in filtered_chain_nodes:

        if node in visited:
            continue

        # Begin of chain
        if left is None :

            chain = []
            current = node

            while current is not None and current not in visited:
                chain.append(current)
                visited.add(current)

                nxt = right_of.get(current)

                current = nxt
            
            if len(chain) > 1 :
                chains.append(chain)

    return chains


def compress_non_branching_paths(graph : Graph):
    """
    Compress chains of non-branching nodes

    Parameter
    ----------
    graph : Graph
        Graph to compress
    """
    
    # Compute non-branching nodes 
    chain_nodes = get_non_branching_nodes(graph)

    # Build non-branching nodes chains
    chains = build_non_branching_paths(chain_nodes)

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
        final_node = chain[0]
        # It will receive a new sequence replacing all nodes of the chain
        new_seq_parts = []
        
        # List for the path having their first node in this chain
        offset_paths = []

        # Course of the chain
        for i in range(len(chain)):

            node = chain[i]

            # If the path begins in the middle of this chain, we treat it later
            if node in first_node_to_paths:
                offset_paths = first_node_to_paths[node]

            # Adding current sequence to new sequence
            new_seq_parts.append(graph.segments[node]['seq'])

            # Nodes of the chain, except first one are supressed
            if i != 0 :
                suppressed_nodes.add(node)
        
        # Build new sequence
        new_seq = ''.join(new_seq_parts)

        # In cases of offset in the chain
        for path_name in offset_paths:
            len_of_chain_in_path = 0
            final_node_in_hap = False
            # Comput which portion of the chain is in the path
            for n, o in graph.paths[path_name]['path']:
                if n not in chain:
                    break
                else :
                    len_of_chain_in_path += graph.segments[n]['length']
                    if n == final_node :
                        final_node_in_hap = True
            # Offset value is the difference between chain length and portion in the path
            len_offset = len(new_seq) - len_of_chain_in_path

            if len_offset != 0 :
                offsets[path_name] = len_offset
            # In cases of offset where the first node of the chain is not in the path
            # The first node of the chain is saved as first node of the path (since onmy the first node of the chain is saved in final graph)
            # Orientation is the one of the real first node of the path
            if not final_node_in_hap:
                first_nodes[path_name] = (final_node, graph.paths[path_name]['path'][0][1])
        
        # Replace content of final node
        graph.segments[final_node]['seq'] = new_seq
        graph.segments[final_node]['length'] = len(new_seq)

    # Compute new nodes set
    new_segments = {}

    for n, data in graph.segments.items():

        if n in suppressed_nodes :
            continue
        
        new_segments[n] = data

    graph.segments = new_segments
    
    # Paths rewriting : only keep the nodes that are still in the new nodes set
    new_paths = {}

    for name, pdata in graph.paths.items():

        new_path = []

        if name in first_nodes:
            new_path.append(first_nodes[name])

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

    for pdata in graph.paths.values():
        for (u, ou), (v, ov) in pairwise(pdata["path"]):
            graph.add_edge(u, ou, v, ov)

    for data in graph.segments.values():
        data["successors"] = []
        data["predecessors"] = []

    # Update
    graph.sequence_offsets()
    graph.compute_neighbors()

    return offsets
