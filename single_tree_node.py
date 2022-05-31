import numpy as np
from tree_functions import loop_children


class SingleTreeNode():

    def __init__(self, x0, y0, width, height, children=None, parents=None, points=None, level=0):

        self.x0 = x0
        self.y0 = y0
        self.w = width
        self.h = height
        self.verts = ((x0, x0 + width), (y0, y0 + height))
        self.center = (x0 + width/2, y0 + height/2)
        self.children = children
        self.cardinal_neighbors = [None, None, None, None,]
        self.nearest_neighbors = None
        self.cardinal_index = 0
        self.parents = parents
        self.points = []
        self.level = level
        self.inner = None 
        self.outer = None

        if points is not None:
            self.add_points(points)

    def __iter__(self):
        if self.has_children():
            for child in self.children:
                yield child

    def __len__(self):
        if self.points is None:
            return 0 
        return len(self.points)

    def add_points(self, points):
        if self.has_children():
            for child in self.children:
                child.add_points(points)
        else:
            for d in points:
                if self.contains(d.x, d.y):
                    self.points.append(d)

    def get_points(self):
        return self.points

    def has_children(self):
        if self.children is not None:
            return self.children

    def get_child(self, i):
        if self.children is None:
            return self
        return self.children[i]

    def split(self):
        if self.has_children():
            return

        w = self.w / 2
        h = self.h / 2
        x0 = self.verts[0][0] 
        y0 = self.verts[1][0]

        # [NW, NE, SW, SE] =  [0, 1, 2, 3]
        self.children = [SingleTreeNode(w, h, xi, yi, points=self.points, level=self.level+1, parent=self)
                        for yi in (y0 + h, y0) for xi in (x0, x0 + w)]

        for i, c in enumerate(self.children):
            c.cardinal_index = i

    def threshold_split(self, threshold):
        if len(self) > threshold:
            self.split()
        if self.has_children():
            for child in self.children:
                child.threshold_split(threshold)

    def contains(self, x, y):
        return ((x >= self.verts[0][0] and x < self.verts[0][1]) and
                (y >= self.verts[1][0] and y < self.verts[1][1]))

    def is_leaf(self):
        if self.children is None:
            return True
        else: 
            return False

    def set_cardinal_neighbors(self):
        for i, child in enumerate(self.children):
            n_sibling = (abs(1 + (i^1) - i), abs(1 + (i^2) - i))

            child.cardinal_neighbors[n_sibling[0]] = self.children[i^1]
            child.cardinal_neighbors[n_sibling[1]] = self.children[i^2]

            n_parent = tuple(set((0,1,2,3)) - set((n_sibling)))
            n_cardinal = lambda j, k: j^((k + 1)%2 + 1)

            child.cardinal_neighbors[n_parent[0]] = (self.cardinal_neighbors[n_parent[0]].get_child(n_cardinal(i, n_parent[1]))
                                        if self.cardinal_neighbors[n_parent[0]] is not None
                                        else None)
            child.cardinal_neighbors[n_parent[1]] = (self.cardinal_neighbors[n_parent[1]].get_child(n_cardinal(i, n_parent[0]))
                                        if self.cardinal_neighbors[n_parent[1]] is not None
                                        else None)
            if child.has_children():
                child.set_cardinal_neighbors()

    def traverse(self):
        if self.has_children():
            for child in loop_children(self):
                yield child

    avoid_index = (1, 2, 0, 3)
    corner_index = (3, 2, 0, 1)

    @property
    def get_nearest_neighbors(self):
        if self.nearest_neighbors is not None:
            return self.nearest_neighbors

        nn = [cn.cardinal_neighbors[(i+1)%4]
              for i, cn in enumerate(self.cardinal_neighbors)
              if cn is not None 
              and cn.level == self.level]

        nn += [cn.cardinal_neighbors[(i+1)%4].get_child(self.corner_index[i])
               for i, cn in enumerate(self.cardinal_neighbors)
               if cn is not None 
               and cn.cardinal_neighbors[(i+1)%4] is not None 
               and (cn.level < self.level and i != self.avoid_index[self.cardinal_index])]

        nn = [n for n in self.cardinal_neighbors + nn if n is not None]
        self.nearest_neighbors = nn
        
        return nn

    def interaction_set(self):
        nn, pn = self.get_nearest_neighbors, self.parent.get_nearest_neighbors
        int_set = []
        for n in pn:
            if n.has_children():
                int_set += [c for c in n if c not in nn]
            elif n not in nn:
                int_set.append(n)
        return int_set