import numpy as np


def distance(point_1, point_2):
    distance = np.sqrt((point_1[0] - point_2[0])**2 + (point_1[1] - point_2[1])**2)
    return distance