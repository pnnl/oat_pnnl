import copy
import numpy as np


def lower_star_adjacency_matrix(dismat=[], densityvec=[], threshold=0):
    """
    # Arguments

    :param dismat: a square symmetric matrix whose entries represent dissimilarity scores
    :param densityvec List[float]: a list of real numbers representing a real value function on vertices
    :param threshold float: a real number

    # Returns

    :return adj: A weighted adjacency matrix where
    - two points are adjacent iff the distance between them is <= `threshold`
    - the weight of an edge connecting two points is the maximum of the values of `densityvec` taken on those two points
    """
    plaid = np.maximum( densityvec.reshape(-1,1), densityvec.reshape(1,-1) )

    deformed = copy.deepcopy(dismat)
    deformed[dismat > threshold] = np.inf
    deformed[ dismat <= threshold ] = plaid[ dismat <= threshold ]

    return deformed   