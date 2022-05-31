import numpy as np


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
        self.card_index = 0
        self.parents = parents
        self.points = []
        self.level = level
        self.inner = None 
        self.outer = None