class Bubble:
    def __init__(self, boundaries):
        self.boundaries = boundaries
        self.parent = None
        self.interior_nodes = set()
        self.size = 0

