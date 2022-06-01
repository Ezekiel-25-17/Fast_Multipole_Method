import numpy as np
from single_tree_node import SingleTreeNode

class QuadTree():

    def __init__(self, points, threshold, bbox=(1,1), boundary='wall'):
        self.threshold = threshold
        self.root = SingleTreeNode(*bbox, 0, 0)
        if boundary == 'periodic':
            self.root.cardinal_neighbors = [self.root, self.root, self.root, self.root]
        elif boundary == 'wall':
            self.root.cardinal_neighbors = [None, None, None, None]
        else:
            raise AttributeError(f'Boundary of type {boundary} not recognized')
        self.build_tree(points)
        self.depth = None

    def __len__(self):
        l = len(self.root)
        for node in self.root.traverse():
            l += len(node)
        return l

    def __iter__(self):
        for points in self.root.get_points():
            yield points

    def build_tree(self, points):
        self.root.add_points(points)
        self.root.threshold_split(self.threshold)
        self.root.set_cardinal_neighbors()

    @property
    def depth(self):
        if self.depth is None:
            self.depth = max([node.level for node in self.root.traverse()])
        return self.depth

    @property
    def nodes(self):
        return [node for node in self.root.traverse()]

    def traverse_nodes(self):
        for node in self.root.traverse():
            yield node


class Point2D():

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.pos = (x, y)


class Particle(Point2D):

    def __init__(self, charge, x, y):
        super(Particle, self).__init__(x, y)
        self.q = charge
        self.phi = 0