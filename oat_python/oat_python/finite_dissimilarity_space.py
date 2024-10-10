import numpy as np
import sklearn
import sklearn.metrics
import scipy

def farthest_point( dismat = np.zeros((0,0)), epsilon = 0 ):
    """
    Use farthest point sampling to select an epsilon net of a disimilarity space.
    The net also contains an epsilon net for every delta > epsilon.

    Returns a sequence of intebers p0, p1, ...
    and a sequence of values e0, e1, ... such that [p0, .., pm]
    is an em net.

    :param dismat: a sparse or dense distance matrix; all off-diagonal zero entries will be treated as infinity
    :param epsilon: a float
    """
    print("This method of farthest point sampling may be slow.  For large point clouds, `try farthest_point_sans_dismat`, which has the added advantage that you don't have to compute a distance matrix.")

    if dismat.shape[0] == 0:
        return np.zeros((0,0))
    if dismat.shape[0] == 1:
        return [0], 0

    def get_disvec(dismat,i):
        """
        If dismat is dense the return dismat[i].
        If dismat is sparse then convert dismat[i] to dense a dense vector disvec
        and replace disvec[j] with infinity whenever disvec[j]=0 (except for j=i)
        """
        if type(dismat) == scipy.sparse._csr.csr_matrix:
            disvec = np.asarray(dismat[i].todense()).reshape(-1)
            disvec[disvec==0] = np.inf
            disvec[i] = 0
            return disvec 
        else:
            return dismat[i]

    disvec          =   get_disvec(dismat,0)  # entry j of this vector represents the min distance from point j to any other point in the point cloud
    maxpnt          =   disvec.argmax()
    maxdis          =   disvec[ maxpnt ]
    maxdiscurve     =   [maxdis] # maxdiscurve[p] = max distance from any point to the first p+1 points of the permutation
    greedyperm      =   [0, maxpnt]      

    while maxdis > epsilon:
        disvec      =   np.minimum( disvec, get_disvec(dismat, maxpnt)  ) # smallest distance to any point in our net
        maxpnt      =   np.argmax( disvec )
        maxdis      =   disvec[ maxpnt ]
        greedyperm.append(maxpnt)
        maxdiscurve.append(maxdis)

    return greedyperm, maxdiscurve

def farthest_point_sans_dismat_fixed_cardinality( cloud, netcardinality=1 ):
    """
    Use farthest point sampling to select an epsilon net of a disimilarity space.
    The net will contain exactly `netcardinality` points.
    """
    if cloud.shape[0] == 0:
        return np.zeros((0,0))
    if cloud.shape[0] == 1:
        return [0], 0

    if netcardinality > cloud.shape[0]:
        print("Error: net cardinality cannot exceed the number of points")
        return 


    # entry j of this vector represents the min distance from point j to any other point in the point cloud    
    disvec      =   sklearn.metrics.pairwise_distances( cloud[[0],:], cloud )
    maxpnt      =   np.argmax( disvec )
    maxdis      =   disvec[ maxpnt ]
    maxdiscurve =   [maxdis] # maxdiscurve[p] = max distance from any point to the first p+1 points of the permutation
    greedyperm  =   [0, maxpnt]    

    for _ in range( netcardinality-1 ):
        disvec_marg =   sklearn.metrics.pairwise_distances( cloud[[maxpnt],:], cloud )
        disvec      =   np.minimum( disvec, disvec_marg  ) # smallest distance to any point in our net
        maxpnt      =   np.argmax( disvec )
        maxdis      =   disvec[ maxpnt ]
        greedyperm.append(maxpnt)
        maxdiscurve.append(maxdis)

    return greedyperm, maxdiscurve 


def farthest_point_sans_dismat( cloud, epsilon=0 ):
    """
    Use farthest point sampling to select an epsilon net of a disimilarity space.
    The net also contains an epsilon net for every delta > epsilon.

    Returns a sequence of intebers p0, p1, ...
    and a sequence of values e0, e1, ... such that [p0, .., pm]
    is an em net.

    :param dismat: a sparse or dense distance matrix; all off-diagonal zero entries will be treated as infinity
    :param epsilon: a float
    """
    if cloud.shape[0] == 0:
        return np.zeros((0,0))
    if cloud.shape[0] == 1:
        return [0], 0
    
    # entry j of this vector represents the min distance from point j to any other point in the point cloud    
    disvec      =   sklearn.metrics.pairwise_distances( cloud[[0],:], cloud )[0]
    maxpnt      =   np.argmax( disvec )
    maxdis      =   disvec[ maxpnt ]      
    maxdiscurve =   [maxdis] # maxdiscurve[p] = max distance from any point to the first p+1 points of the permutation
    greedyperm  =   [0, maxpnt]    

    while maxdis > epsilon:
        disvec_marg =   sklearn.metrics.pairwise_distances( cloud[[maxpnt],:], cloud )[0]
        disvec      =   np.minimum( disvec, disvec_marg  ) # smallest distance to any point in our net
        maxpnt      =   np.argmax( disvec )
        maxdis      =   disvec[ maxpnt ]
        greedyperm.append(maxpnt)
        maxdiscurve.append(maxdis)

    return greedyperm, maxdiscurve    


def hausdorff_distance( dissimilarity_matrix=None ):
    """
    Returns max( min( dissimilarity_matrix(x_j, X_i) ) ) for x_j in X_j,  where i \neq j
    """
    # pdist = sklearn.metrics.pairwise_distances(pc1,pc2)
    a = np.max( np.min( dissimilarity_matrix, axis=0 ) )
    b = np.max( np.min( dissimilarity_matrix, axis=1 ) )
    return np.maximum( a, b )


