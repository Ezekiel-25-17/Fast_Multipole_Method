from single_tree_node import SingleTreeNode
import numpy as np


class QuadTree():

    def __init__(self, points, threshold, bbox=(1,1), boundary='wall'):
        self.threshold = threshold
        self.root = SingleTreeNode(*bbox, 0, 0)
        if boundary == 'periodic':
            self.root.cardinal_neighbors = 4*[self.root,]
        elif boundary == 'wall':
            self.root.cardinal_neighbors = 4*[None,]
        else:
            raise AttributeError(f'Boundary of type {boundary} not recognized')
        self.build_tree(points)
        self._depth = None

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
        if self._depth is None:
            self.depth = max([node.level for node in self.root.traverse()])
        return self._depth

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

    def __init__(self, x, y, charge):
        super(Particle, self).__init__(x, y)
        self.q = charge
        self.phi = 0


eps = (7/3 - 4/3) -1

def build_tree(points, tree_threshold=None, bbox=None, boundary='wall'):
    if bbox is None:
        coords = np.array([(p.x, p.y) for p in points])
        bbox = (max(coords[:, 0]) + eps, max(coords[:, 1]) + eps)
    if tree_threshold is None:
        tree_threshold = 5

    return QuadTree(points, tree_threshold, bbox=bbox, boundary=boundary)