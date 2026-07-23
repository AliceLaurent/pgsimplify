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

def filters_nodes_being_path_ends(nodes, end_nodes):
    """
    Filter nodes that are path ends which prevent compression

    Parameter
    ----------
    nodes : List[(str, str, str)]
        The nodes of the form (left, node, right) with left and right its potential neighbours in a chain
    end_nodes : dict[set]
        Gives for each node being at the end of one or more paths the nodes just before/after 
    
    Returns
    -------
    List[(str, str, str)]
        List of filtered nodes
    """
    filtered_nodes = []
    for node in nodes :
        (left, n, right) = node
        if n in end_nodes:
            # The part of the chain on the side where the paths go can still be compressed
            left_possible = True
            right_possible = True
            for snd in end_nodes[n]:
                if left == snd :
                    right_possible = False
                elif right == snd:
                    left_possible = False
                else :
                    left_possible = False
                    right_possible = False
            if left_possible and right_possible :
                filtered_nodes.append(node)
            elif left_possible :
                filtered_nodes.append((left, n, None))
            elif right_possible :
                filtered_nodes.append((None, n, right))
        else :
            filtered_nodes.append(node)
    return filtered_nodes
                

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

    # For cases where a path begins in the middle of a future compressed node
    end_nodes = defaultdict(set)

    # Build dictionnary end node of a path -> node just before/after end node for direction
    for name, pdata in graph.paths.items():
        end_nodes[pdata['path'][0][0]].add(pdata['path'][1][0])
        end_nodes[pdata['path'][-1][0]].add(pdata['path'][-2][0])

    filtered_nodes = filters_nodes_being_path_ends(chain_nodes, end_nodes)

    # Build non-branching nodes chains
    chains = build_non_branching_paths(filtered_nodes)

   # Compute nodes that will be removed and new sequence that will replace them for each chain
    suppressed_nodes = set()
    for chain in chains:

        # First node of the chain will replace the entire chain
        final_node = chain[0]
        # It will receive a new sequence replacing all nodes of the chain
        new_seq_parts = []

        # Course of the chain
        for i in range(len(chain)):

            node = chain[i]


            # Adding current sequence to new sequence
            new_seq_parts.append(graph.segments[node]['seq'])

            # Nodes of the chain, except first one are supressed
            if i != 0 :
                suppressed_nodes.add(node)
        
        # Build new sequence
        new_seq = ''.join(new_seq_parts)
        
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
