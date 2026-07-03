"""
Structure to represent the snarl structure.
Unlike usual representation containing both snarl nodes and snarl chains boundary nodes, there is only one type of node representing snarls
With this structure, the snarl structure is a snarl forest, because nothing holds first level snarls together (no snarl chain represented)
An artificial root node is added to link all trees, we can imagine the whole graph is included in one root snarl 
"""
import json


class SnarlNode:
    def __init__(self, key, parent_key=None):
        self.key = key
        self.parent_key = parent_key
        self.parent = None
        self.children = []


class SnarlTree:
    def __init__(self, json_path: str):
        self.json_path = json_path
        self.nodes = {}
        # Artificial root node
        self.root = SnarlNode(key=(-1, -1), parent_key=None)

        self.index_to_key = {}
        self.key_to_index = {}
        self._parse_json()
        self._build_tree()
        self.flatten()

    # Sort pairs to guarantee unicity
    def _canonical_pair(self, a, b):
        a, b = sorted([int(a), int(b)])
        return (a, b)

    # Parse json snarl file
    def _parse_json(self):

        with open(self.json_path) as f:
            for line in f:

                line = line.strip()
                if not line:
                    continue

                data = json.loads(line)

                start = data["start"]["node_id"]
                end = data["end"]["node_id"]

                key = self._canonical_pair(start, end)

                parent_key = None
                if "parent" in data:
                    p = data["parent"]
                    parent_key = self._canonical_pair(
                        p["start"]["node_id"],
                        p["end"]["node_id"]
                    )

                self.nodes[key] = SnarlNode(key, parent_key)

    # Build the snarl tree using dictionnary build during json parsing
    def _build_tree(self):

        waiting = list(self.nodes.values())
        prev_size = -1

        while waiting and len(waiting) != prev_size:

            prev_size = len(waiting)
            new_waiting = []

            for node in waiting:

                # pas de parent → ROOT
                if node.parent_key is None:
                    node.parent = self.root
                    self.root.children.append(node)
                    continue

                parent = self.nodes.get(node.parent_key)

                if parent and parent.parent is not None:
                    node.parent = parent
                    parent.children.append(node)
                else:
                    new_waiting.append(node)

            waiting = new_waiting


    # praitn snarl tree in a text file
    def print_tree(self):

        def rec(node, prefix="", last=True):

            connector = "└── " if last else "├── "

            if node.key == (-1, -1):
                print("ROOT")
            else:
                print(prefix + connector + f"{node.key[0]}--{node.key[1]}")

            prefix += "    " if last else "│   "

            for i, child in enumerate(node.children):
                rec(child, prefix, i == len(node.children) - 1)

        rec(self.root)

    # Flatten the tree to obtain all snarls and compute dictionnary to access snarl keys and ans index
    def flatten(self):

        index_to_key = {}
        key_to_index = {}

        stack = [self.root]
        idx = 0

        while stack:

            node = stack.pop()

            # Ignoring artificial root
            if node.key != (-1, -1):
                index_to_key[idx] = node.key
                key_to_index[node.key] = idx
                idx += 1

            for child in reversed(node.children):
                stack.append(child)

        self.index_to_key = index_to_key
        self.key_to_index = key_to_index