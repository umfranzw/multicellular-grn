class DrawTree():
    def __init__(self, cell, data_tools, parent=None, depth=0, number=1):
        self.x = -1.
        self.y = depth
        self.cell = cell
        
        self.children = []
        for i, child in enumerate(data_tools.get_cell_children(cell)):
            self.children.append(DrawTree(child, self, depth + 1, i + 1))
            
        self.parent = parent
        self.thread = None
        self.mod = 0
        self.ancestor = self
        self.change = self.shift = 0
        self._lmost_sibling = None
        #this is the number of the node in its group of siblings 1..n
        self.number = number

    def left(self): 
        return self.thread or len(self.children) and self.children[0]

    def right(self):
        return self.thread or len(self.children) and self.children[-1]

    def lbrother(self):
        n = None
        if self.parent:
            for node in self.parent.children:
                if node == self: return n
                else:            n = node
        return n

    def get_lmost_sibling(self):
        if not self._lmost_sibling and self.parent and self != \
        self.parent.children[0]:
            self._lmost_sibling = self.parent.children[0]
        return self._lmost_sibling
    
    lmost_sibling = property(get_lmost_sibling)

    def __str__(self): return "%s: x=%s mod=%s" % (self.cell, self.x, self.mod)
    def __repr__(self): return self.__str__()
