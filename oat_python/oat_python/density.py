import copy
import numpy as np 

import sklearn
from sklearn.neighbors import NearestNeighbors

def density( dissimilarity_matrix ):
    """
    :param dissimilarity_matrix: a square matrix representing pairwise dissimilarity between data points
    :return densities: a numpy array representing the "density" score of each data point
    """
    densities = []
    dissimilarity_matrix = dissimilarity_matrix / np.mean( dissimilarity_matrix )
    for vec in dissimilarity_matrix:
        densities.append( np.mean( (1+vec)**(-8 ) ) )
    return np.array(densities)    

def density_nn( cloud, nn=100 ):
    """
    Compute density for each point according to its `nn` nearest neighbors.

    First we normalize the point cloud by variance, the compute density.

    :param cloud: a point cloud
    :param nn (default 100): number of neighbors to use in the density calculation
    :return densities: a numpy array representing the "density" score of each data point
    """

    # calculate variance (for normalization purposes)
    var = np.var( cloud, axis=0 ).sum()

    # nbrs = NearestNeighbors(n_neighbors=100, algorithm='ball_tree').fit(cloudtest)
    nbrs = NearestNeighbors(n_neighbors=nn ).fit(cloud)    
    distances, indices = nbrs.kneighbors(cloud)
    return np.sum( (1+distances/var)**(-8), axis=1 ).reshape(-1)
