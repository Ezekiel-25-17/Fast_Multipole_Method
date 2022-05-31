import numpy as np
from quadtree import QuadTree

def loop_children(parent):
    for child in parent.children:
        if child.children:
            for subchild in loop_children(child):
                yield subchild
        yield child


def build_tree(points, tree_threshold=None, eps = (7/3)-(4/3)-1, bbox=None, boundary='wall'):
    if bbox is None:
        coords = np.array([(p.x, p.y) for p in points])
        bbox = (max(coords[:, 0]) + eps, max(coords[:, 1]) + eps)
    if tree_threshold is None:
        tree_threshold = 5

    return QuadTree(points, tree_threshold, bbox=bbox, boundary=boundary)
