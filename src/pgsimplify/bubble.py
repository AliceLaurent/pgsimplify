class Bubble:
    def __init__(self, boundaries):
        self.boundaries = boundaries
        self.parent = None
        self.interior_nodes = set()
        self.paths = []
        self.size = 0

