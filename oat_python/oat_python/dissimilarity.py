"""
A module for dissimilarity data
"""

#   SEE THE BOTTOM OF THIS FILE FOR UNIT TESTS

from cmath import inf
import numpy as np
import sklearn
import sklearn.metrics
from . import point_cloud
import scipy
from sklearn.neighbors import radius_neighbors_graph
import itertools
import warnings

# Ignore a specific warning
warnings.filterwarnings("ignore", category=scipy.sparse.SparseEfficiencyWarning)




# class FormattedDissimilarityMatrix:
#     def __init__(self, sparse_matrix):
#         """
#         Wraps a `sparse_matrix` in a `FormattedDissimilarityMatrix`.

#         The constructor will check that the input is a symmetric (with zero margin of error) scipy sparse `csr_matrix`.

#         It will also store explicit an explicit zero value in entry `[p,p]`, for all `p` such that `sparse_matrix[p,p]==0`.
#         """

#         print("\n\nDirect construction of FormattedDissimilarityMatrix often yields suboptimal results; consider the other constructors offered in the `dissimilarity` module. \n\n")

#         if not isinstance(sparse_matrix, scipy.sparse.csr_matrix):
#             raise TypeError("Input matrix must be a scipy.sparse.csr_matrix.")
#         transposed = sparse_matrix.transpose
#         if not all( [transposed.data == sparse_matrix.data, transposed.indptr == sparse_matrix.indptr, transposed.indices == sparse_matrix.indices ] ):
#             raise Exception("Input matrix must be symmetric")  
              
#         # ensure that all diagonal entries are stored explicitly
#         missing_diagonal_indices = [p for p in range(sparse_matrix.shape[0]) if sparse_matrix[p,p] == 0 ]
#         sparse_matrix[ missing_diagonal_indices, missing_diagonal_indices ] = 0

#         self._dissimilarity_matrix = sparse_matrix

#     def get_matrix(self):
#         """
#         Returns the internally stored sparse matrix
#         """
#         return self._distance_matrix        

def enclosing_from_cloud_slow( cloud, argminmax=False ):
    """
    Calculates the enclosing radius of a point cloud without calculating a dense distance matrix.

    This method can be paired with `oat.dissimilarity.matrix_from_cloud_slow` in order to generate a
    sparse dissimilarity matrix with the same Vietoris Rips persistent homology as the dense
    dissimilarity matrix as the point cloud.

    # Definition

    The enclosing radius of a point cloud is obtained from its distance matrix `D` by (i) taking the maximum
    of each row of `D`, then (ii) taking the minima of these values, over all rows.

    The enclosing radius produced by this matrix is **guaranteed to be compatible** with the distance matrix
    produced by `matrix_from_cloud_slow`, in the sense that `matrix_from_cloud_slow( cloud, dissimilarity_masx = enclosing_from_cloud_slow(cloud) )`
    has the same Vietoris Rips persistent homology as `matrix_from_cloud_slow( cloud, dissimilarity_masx = inf )`

    # Arguments

    :param cloud array-like: a point cloud represented as a list of points, e.g. a
    list of tuples or a numpy array of size num_points x dimension; each slice `cloud[i]`
    will be treated as a point
    
    # Returns

    :return enclosing radius float: if `argminmax=False`, returns the enclosuing radius
    :return enclosing_radius dict: if `argminmax=True`, returns  a dictionary with enclosing radius and the indices of the two points where it occurs
    """

    # if the cloud is empty, then we're taking a minimum over an empty set, which is inf
    if np.shape(cloud)[0] == 0:
        return inf
    
    # if the cloud lives in R^0 then every point is distance 0 from every other point
    if np.shape(cloud)[1] == 0:
        return 0
    

    def argmaxima_iterator():
        """
        For each point, yields the index of the most distant point
        """
        for point in cloud:
            yield max( range(len(cloud)), key= lambda i: euclidean_distance( cloud[i], point ) )

    ( row, col ) = min( enumerate( argmaxima_iterator() ), key = lambda x : euclidean_distance( cloud[x[0]], cloud[x[1]] ) )

    enclosing_radius = euclidean_distance(cloud[row], cloud[col])

    if argminmax:
        return dict( pointa=row, pointb=col, enclosing_radius=euclidean_distance(cloud[row], cloud[col]))
    else:
        return enclosing_radius
    

def enclosing_from_cloud( cloud, argminmax=False ):
    """
    Calculates the enclosing radius of a point cloud without calculating a dense distance matrix.

    Subject to numerical error; to compensate, one will have to add a small amount to the output, e.g. 0.000000001.
    """

    # if the cloud is empty, then we're taking a minimum over an empty set, which is inf
    if np.shape(cloud)[0] == 0:
        return inf
    
    # if the cloud lives in R^0 then every point is distance 0 from every other point
    if np.shape(cloud)[1] == 0:
        return 0

    npoints     =   np.shape(cloud)[0]
    tree        =   sklearn.neighbors.KDTree(cloud)
    minmax_val  =   np.inf
    minmax_row  =   0
    for row in range(cloud.shape[0]):
        (cols,vals) =   tree.query_radius(  
                            np.array(cloud[row]).reshape(1,-1), 
                            minmax_val + 0.0000000001, 
                            return_distance=True 
                        )
        cols    =   cols[0] # index into first (and only) row of this 2d numpy array
        vals    =   vals[0] # index into first (and only) row of this 2d numpy array
        if len(cols) == npoints: 
            col_rel         =   vals.argmax() 
            val             =   vals[col_rel]
            if val < minmax_val:
                minmax_val  =   val
                minmax_row  =   row
                minmax_col  =   cols[col_rel]

    if argminmax:
        cols        =   tree.query_radius(  
                                np.array(cloud[minmax_row]).reshape(1,-1), 
                                minmax_val + 0.0000000001, 
                                return_distance=True,  
                        )[-1] 
        max_ind     =   cols[1][0].argmax()
        minmax_col  =   cols[0][0][max_ind]
        # and take the last column in the sequence
     

        return dict( pointa=minmax_row, pointb=minmax_col, enclosing_radius=minmax_val)
    else:
        return minmax_val    
    

def enclosing_from_csr( matrix, argminmax=False ):
    """
    The enclosing radius of a dissimilarity matrix, quantified as the minimum of the maxima of the rows.

    - Entries without explicit stored values are treated as `inf`.
    - If the matrix has no rows, then `inf` is returned.
    - If the matrix has at least one row but no columns, then `-inf` is returned.

    # Definition

    The enclosing radiuos of (distance) matrix `min(v)`, where `v[i]` is the maximum of the values in row `i`.
    If the matrix has no rows, then `inf` is returned.

    # Significance

    The enclosing radius is significant for Vietoris Rips persistence calculations, because the homology
    of a filtered Vietoris Rips complex is zero for filtration parameter values equal to and exceeding the
    enclosing radius.

    # Caution (numerical error and asymmetry)

    Distance matrices produces in Python often suffer from numerical error; for example,

    - `sklearn.neighbors.radius_neighbors_graph` often produces asymmetric matrices
    - `sklearn.metrics.pairwise_distances` often produces asymmetric matrices

    Consequently, calculating enclosing radius with this function, and using it as a threshold for
    excluding entries in a dissimilarity matrix, may yield a sparse dissimilarity matrix with the
    wrong persistent homology.  This can be avoided by adding a tiny amoung of buffer to the
    enclosuing radius, e.g. `enclosing_radius + 0.0000000001`.
    
    # Arguments

    :param matrix scipy.sparse.csr
    
    # Returns

    :return enclosing float: if `argminmax = False`, returns the enclosing radius

    :return enclosing_radius dict: if `argminmax = True`, returns a dictionary with enclosing radius and the row/column where it occurs
    """
    m,n = matrix.shape

    if m == 0:
        if argminmax:
            return dict( row=None, col=None, enclosing_radius=inf)
        else:
            return inf
    elif n == 0:
        if argminmax:
            return dict( row=None, col=None, enclosing_radius=-inf)
        else:
            return -inf

    def argmaxima_iterator():
        """
        For each row, yields the index of the column where the maximum value occurs
        """
        for row_index in range(m):
            linear_indices =  range(matrix.indptr[row_index], matrix.indptr[row_index+1])
            if len(linear_indices) < n:
                yield inf
            else:
                index_max = max( linear_indices, key= lambda i: matrix.data[i] )
                yield matrix.data[index_max]

    (row, enclosing_radius ) = min( enumerate(argmaxima_iterator()), key= lambda x: x[1] )

    col = matrix[row].argmax()

    if argminmax:
        return dict( row=row, col=col, enclosing_radius=enclosing_radius )
    else:
        return enclosing_radius
    


def enclosing_from_dense( matrix, argminmax=False ):
    """
    The enclosing radius of a dense dissimilarity matrix, quantified as the minimum of the maxima of the rows.

    - If the matrix has no rows, then `inf` is returned.
    - If the matrix has at least one row but no columns, then `-inf` is returned.

    # Definition

    The enclosing radiuos of (distance) matrix `min(v)`, where `v[i]` is the maximum of the values in row `i`.

    # Significance

    The enclosing radius is significant for Vietoris Rips persistence calculations, because the homology
    of a filtered Vietoris Rips complex is zero for filtration parameter values equal to and exceeding the
    enclosing radius.

    # Caution (numerical error and asymmetry)

    Distance matrices produces in Python often suffer from numerical error; for example,

    - `sklearn.neighbors.radius_neighbors_graph` often produces asymmetric matrices
    - `sklearn.metrics.pairwise_distances` often produces asymmetric matrices

    Consequently, calculating enclosing radius with this function, and using it as a threshold for
    excluding entries in a dissimilarity matrix, may yield a sparse dissimilarity matrix with the
    wrong persistent homology.  This can be avoided by adding a tiny amoung of buffer to the
    enclosuing radius, e.g. `enclosing_radius + 0.0000000001`.
    
    # Arguments

    :param matrix scipy.sparse.csr
    
    # Returns

    :return enclosing float: if `argminmax = False`, returns the enclosing radius

    :return enclosing_radius dict: if `argminmax = True`, returns a dictionary with enclosing radius and the row/column where it occurs
    """
    m,n = matrix.shape

    if m == 0:
        if argminmax:
            return dict( row=None, col=None, enclosing_radius=inf)
        else:
            return inf
    elif n == 0:
        if argminmax:
            return dict( row=None, col=None, enclosing_radius=-inf)
        else:
            return -inf

    row     =   matrix.max(axis=1).argmin()
    col     =   matrix[row].argmax()
    enclosing_radius = matrix[row][col]

    if argminmax:
        return dict( row=row, col=col, enclosing_radius=enclosing_radius )
    else:
        return enclosing_radius    



def matrix_from_cloud( cloud, dissimilarity_max ):
    """
    Returns a sparse (Scipy CSR) Euclidean distance matrix, where all entries with value strictly greater than `dissimilarity_max` are dropped.

    # Arguments

    :param cloud: any format for a point cloud compatible with `sklearn.neighbors.radius_neighbors_graph`
    :param dissimilarity_max float: a non-negative real number; all distances with value above this threshold are dropped

    # Assymetry

    This function uses `sklearn.neighbors.radius_neighbors_graph` to construct a dissimilarity matrix, however
    the output of this function is **not symmetric** in general.  Therefore the construction process takes two steps:

    - construct a matrix `A` with `sklearn.neighbors.radius_neighbors_graph`
    - replace `A` with the entrywise-maximum of `A` and the transpose of `A`
    - store explicit zero values for all diagonal entries
    - if `dissimilarity_max` is nonnegative, then store explicit zero values for all diagonal entries;
      otherwise return an empty sparse matrix of appropriate size     

    # Examples

    ```
    import numpy as np
    import oat_python as oat

    cloud           =   np.random.rand(10,2)
    dissimilarity   =   oat.dissimilarity.matrix_from_cloud_slow(
                            cloud               =   cloud, 
                            dissimilarity_max   =   oat.dissimilarity.enclosing_radius(cloud) 
                        )
    ``` 

    """

    # if the cloud is empty then return a 0 x 0 matrix
    if np.shape(cloud)[0] == 0:
        return scipy.sparse.csr_matrix((0,0))

    if dissimilarity_max < 0:
        num_rows = np.shape(cloud)[0] # number of points in the cloud
        return scipy.sparse.csr_matrix((num_rows,num_rows)) # return an empty matrix of appropriate size

    A   =   radius_neighbors_graph( 
                cloud, 
                radius = dissimilarity_max, 
                mode='distance', 
                include_self=True
            )    
    A   =   A.maximum( A.T ) # ensure the matrix is symmetric
    
    A.setdiag(0) # set all diagonal entries to zero; otherwise the matrix will be empty, including along the diagonal, which is correct behavior
    A.sort_indices()
    
    return A

#   DEPRECATED; OK TO DELETE
# def sparse_distance_matrix(
#         cloud,
#         distance_max            =   None,
#         enclosing_tolerance     =   None,
#     ):
#     """
#     Returns a sparse, symmetric`scipy.sparse.csr` matrix representing
#     the pairwise distance between points.  Distances over the provided
#     threshold are excluded.

#     If `distance_max` is `None`, then the enclosing radius of the point
#     cloud will be used.
#     """

#     if distance_max is None:
#         if enclosing_tolerance is None:
#             A           =   scipy.sparse.csr_matrix( sklearn.metrics.pairwise_distances(cloud) )
#             enclosing   =   None
#         else:            
#             enclosing   =   enclosing_radius(cloud)
#             A           =   radius_neighbors_graph( 
#                                 cloud, 
#                                 enclosing["enclosing_radius"] * (1+enclosing_tolerance), 
#                                 mode='distance', 
#                                 include_self=True
#                             )  
#     else:
#         if not enclosing_tolerance is None:
#             Exception("Either `eclosing_tolerance` or `distance_max` must be `None`")
#         else:
#             A           =   radius_neighbors_graph( 
#                                 cloud, 
#                                 distance_max, 
#                                 mode='distance', 
#                                 include_self=True
#                             )       
#             enclosing   =   None    
#     A   =   A.maximum( A.transpose() ) # ensures the matrix is symmetric
#     A[ range(len(cloud)), range(len(cloud)),  ] = 0
#     A.sort_indices()

#     return dict(distance_matrix = A, enclosing=enclosing)


def matrix_from_dense( dissimilarity_matrix, dissimilarity_max = np.inf ):
    """
    Converts a dense dissimilarity matrix to a scipy sparse CSR matrix that meets
    the formatting requirements of the OAT persistent homology solver.
    
    In particular, this constructor checks that

    - the input is a
      - square
      - symmetric
      - numpy array
    - for each i, the entry (i,i) takes the smallest value of any entry in row i
    
    It then stores the input in a scipy sparse CSR matrix, where
      - all entries with value strictly greater than `dissimilarity_max` are removed
      - all entries with value less than or equal to `dissimilarity_max` are stored explicitly -- including entries with value 0
    """

    if not isinstance( dissimilarity_matrix, np.ndarray ):
        raise TypeError("The input to `oat.dissimilarity.matrix_from_dense` must be a numpy array")

    if (len(dissimilarity_matrix.shape)!=2) or ( dissimilarity_matrix.shape[0] != dissimilarity_matrix.shape[1]):
        raise Exception(f"Input must be 2-dimensional square matrix, but the matrix provided has shape {dissimilarity_matrix.shape}")
    
    # if the matrix has size 0, then return the size-zero sparse matrix
    if dissimilarity_matrix.shape[0] == 0:
        return scipy.sparse.csr_matrix((0,0))
    
    for p in range(dissimilarity_matrix.shape[0]):
        if dissimilarity_matrix[p,p] != dissimilarity_matrix[p].min():
            raise Exception("The diagonal entry in each row must take the minimum value in that row.")
        for q in range(p, dissimilarity_matrix.shape[0]):
            if not dissimilarity_matrix[p,q] == dissimilarity_matrix[q,p]:
                raise Exception("Input must be symmetric (with zero margin of error)")
    
    data = []; row=[]; col=[];
    for ((i,j),val) in np.ndenumerate( dissimilarity_matrix ):
        if val <= dissimilarity_max:
            data.append(val)
            row.append(i)
            col.append(j)
    
    matrix = scipy.sparse.csr_matrix( 
                ( data, ( row, col ) ), 
                shape = dissimilarity_matrix.shape 
            )
    matrix.sort_indices()
    return matrix


def matrix_from_csr( dissimilarity_matrix, dissimilarity_max = np.inf, set_missing_diagonal_entries_to_zero=True ):
    """
    Converts a scipy sparse CSR dissimilarity matrix to a scipy sparse COO matrix that meets
    the formatting requirements of the OAT persistent homology solver.
    
    In particular, this constructor checks that

    - the input is a
      - square
      - symmetric
      - scipy.sparse.csr_matrix
    - for each i, if row i contains at least one explicit entry then entry (i,i)
      has an explicty entry, and it is the minimum of all explicit entries stored in
      row i
    
    It then stores the input in a scipy sparse COO matrix, where
      - all entries with value strictly greater than `dissimilarity_max` are removed
      - all entries with value less than or equal to `dissimilarity_max` are stored explicitly -- including entries with value 0

    If `set_missing_diagonal_entries_to_zero=True`, then replaces every missing diagonal entry 
    """

    if not isinstance( dissimilarity_matrix, scipy.sparse.csr_matrix ):
        raise TypeError("The input to oat.dissimilarity.matrix_from_csr must be a `scipy.sparse.csr_matrix`")
    
    if dissimilarity_matrix.shape[0] != dissimilarity_matrix.shape[1]:
        raise Exception(f"Input must be square, but the matrix provided has shape {dissimilarity_matrix.shape}")    
    
    for p in range(dissimilarity_matrix.shape[0]):
        row = dissimilarity_matrix.getrow(p)
        if len(row.data) > 0:
            if not p in row.indices:
                raise Exception(f"\n\nRow {p} contains an explicitly stored entry, but no explicitly stored entry in column {p}.\n\nIn most cases users want all diagonal entries to equal zero, so you can fix this problem by calling `matrix.setdiag(0)`.\n\n")
            if dissimilarity_matrix[p,p] > np.min( row.data ):
                raise Exception(f"Entry {(p,p)} contains an explicit entry, but it's not the smallest explicit entry in row {p}")

        for q in range(p, dissimilarity_matrix.shape[0]):
            if not dissimilarity_matrix[p,q] == dissimilarity_matrix[q,p]:
                raise Exception("Input must be symmetric (with zero margin of error)")
    
    dissimilarity_matrix      =   scipy.sparse.coo_matrix(dissimilarity_matrix)
    I           =   [ i for (i, val) in enumerate(dissimilarity_matrix.data) if val <= dissimilarity_max ]
    dissimilarity_matrix.row  =   dissimilarity_matrix.row[I]
    dissimilarity_matrix.col  =   dissimilarity_matrix.col[I]
    dissimilarity_matrix.data =   dissimilarity_matrix.data[I]
    dissimilarity_matrix      =   scipy.sparse.csr_matrix(dissimilarity_matrix)
    dissimilarity_matrix.sort_indices()
    
    return dissimilarity_matrix


def farthest_point( dissimilarity_matrix = np.zeros((0,0)), epsilon = 0 ):
    """
    Use farthest point sampling to select an epsilon net of a disimilarity space.
    The net also contains an epsilon net for every delta > epsilon.

    Returns a sequence of intebers p0, p1, ...
    and a sequence of values e0, e1, ... such that [p0, .., pm]
    is an em net.

    :param dissimilarity_matrix: a sparse or dense distance matrix; all off-diagonal zero entries will be treated as infinity
    :param epsilon: a float
    """
    print("This method of farthest point sampling may be slow.  For large point clouds, `try farthest_point_with_cloud`, which has the added advantage that you don't have to compute a distance matrix.")

    if dissimilarity_matrix.shape[0] == 0:
        return np.zeros((0,0))
    if dissimilarity_matrix.shape[0] == 1:
        return [0], 0

    def get_disvec(dissimilarity_matrix,i):
        """
        If dissimilarity_matrix is dense the return dissimilarity_matrix[i].
        If dissimilarity_matrix is sparse then convert dissimilarity_matrix[i] to dense a dense vector disvec
        and replace disvec[j] with infinity whenever disvec[j]=0 (except for j=i)
        """
        if type(dissimilarity_matrix) == scipy.sparse._csr.csr_matrix:
            disvec = np.asarray(dissimilarity_matrix[i].todense()).reshape(-1)
            disvec[disvec==0] = np.inf
            disvec[i] = 0
            return disvec 
        else:
            return dissimilarity_matrix[i]

    disvec          =   get_disvec(dissimilarity_matrix,0)  # entry j of this vector represents the min distance from point j to any other point in the point cloud
    maxpnt          =   disvec.argmax()
    maxdis          =   disvec[ maxpnt ]
    maxdiscurve     =   [maxdis] # maxdiscurve[p] = max distance from any point to the first p+1 points of the permutation
    greedyperm      =   [0, maxpnt]      

    while maxdis > epsilon:
        disvec      =   np.minimum( disvec, get_disvec(dissimilarity_matrix, maxpnt)  ) # smallest distance to any point in our net
        maxpnt      =   np.argmax( disvec )
        maxdis      =   disvec[ maxpnt ]
        greedyperm.append(maxpnt)
        maxdiscurve.append(maxdis)

    return greedyperm, maxdiscurve

def farthest_point_with_cloud_fixed_cardinality( cloud, netcardinality=1 ):
    """
    Use farthest point sampling to select an epsilon net of a disimilarity space.
    The net will contain exactly `netcardinality` points.
    """
    if cloud.shape[0] == 0:
        return np.zeros((0,0))
    if cloud.shape[0] == 1:
        return [0], 0

    if netcardinality > cloud.shape[0]:
        Exception("Error: net cardinality cannot exceed the number of points")


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


def farthest_point_with_cloud( cloud, epsilon=0 ):
    """
    Use farthest point sampling to select an epsilon net of a disimilarity space.
    The net also contains an epsilon net for every delta > epsilon.

    Returns a sequence of intebers p0, p1, ...
    and a sequence of values e0, e1, ... such that [p0, .., pm]
    is an em net.

    :param dissimilarity_matrix: a sparse or dense distance matrix; all off-diagonal zero entries will be treated as infinity
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






#   ==================================================
#   POINT CLOUD
#   ==================================================





def euclidean_distance_onesided( pointa, pointb ):
    """
    Returns the Euclidean distance between two points, but reversing point order may yield different results due to numerical error
    
    :param pointa array-like
    :param pointb array-like    
    """

    if np.prod(np.shape(pointa)) != np.prod(np.shape(pointb)):
        raise Exception("Error:points must have equal numbers of elements")

    sum_of_square_differences = 0
    for (p,q) in zip( np.nditer(pointa), np.nditer(pointb) ):
        sum_of_square_differences += ( p - q )**2
    
    return np.sqrt( sum_of_square_differences )  

def euclidean_distance( pointa, pointb ):
    """
    Returns the Euclidean distance between two points, 
    """    
    return max( euclidean_distance_onesided(pointa,pointb), euclidean_distance_onesided(pointb,pointa) )

def distance_iterator( cloud, reference_point ):
    """
    Iterates over the Euclidean distances from a point to every point in a point cloud

    :param reference_point array-like
    :param cloud array-like 
    """
    for point in cloud:
        if np.prod(np.shape(reference_point)) != np.prod(np.shape(cloud)[1:]):
            raise Exception("Error: a point in the point cloud does not have the same dimension as the reference point")        
        yield euclidean_distance( point, reference_point )  


#   DEVELOPERS; THIS WAS ORIGINALLY INTENDED TO CIRCUMVENT THE PROBLEM OF NUMERICAL ERROR IN sklearn's radius_neighbors
#               FUNCTION, BUT THAT FUNCTION USES SPECIAL METHODS TO ACCELERATE COMPUTATION, WHICH WE MAY NOT BE ABLE TO
#               MATCH; WE'LL HAVE TO FINISH DEFINING THIS FUNCTION, AND SEE
#   
def matrix_from_cloud_slow( cloud, dissimilarity_max ):
    """
    Returns a symmetric Euclidean distance matrix, where entries strictly above `dissimilarity_max` are not explicitly stored.

    # Performance and numerical error

    This method uses a brute-force function to calculate the distance matrix.  The `sklearn` function `radius_neighbors`
    may be able to give much better performance, and we provide a convenience
    function `matrix_from_cloud`, to format the output in a manner compatible with the OAT persistence
    solver.    

    **However** due to problems of numerical error, the nearest neighbor method won't always produce a matrix equal to the
    one returned by this function.  In particular, calling the nearest neighbors approach with `dissimilarity_max`
    equal to the enclosing radius of the point cloud may yield a dissimilarity matrix whose persistent homology is **not equal**
    to that of the point cloud.  This can be remedied by adding a tiny margin of error to the enclosing radius, e.g.
    `dissimilarity_max = enclosing_radius + 0.00000001`.
    """

    if np.ndim(cloud) ==0:
        Exception("The input to `distance_matrix` must be an array-like object of dimension at least 1")
    
    shape = np.shape( cloud )
    if dissimilarity_max < 0:
        return scipy.sparse.csr_matrix((shape[0],shape[0]))
    
    data    =   []
    row     =   []
    col     =   []
    for row_num in range(shape[0]):
        data.append(0); row.append(row_num); col.append(row_num);
        row_vec = cloud[row_num]
        for col_num in range(row_num+1, shape[0]):
            x = euclidean_distance( row_vec, cloud[col_num] )
            if x <= dissimilarity_max:
                data.append(x)  # add an entry
                row.append(row_num)
                col.append(col_num)
                data.append(x)  # add its transpose
                row.append(col_num)
                col.append(row_num)   
    matrix  =   scipy.sparse.coo_matrix( 
                    (
                        data, ( row, col, ),
                    ),
                    shape = (shape[0],shape[0]),                    
                )
    matrix  =   scipy.sparse.csr_matrix(matrix)
    matrix.sort_indices()
    return matrix
            






#   ===================


def assert_almost_equal_csr( a, b, decimal, err_msg='' ):
    """
    Checks that a and b have the same sparsity pattern and the same explicitly stored entries, up to tolerance.

    Throws an error if the test fails
    """
    np.testing.assert_equal(a.indptr, b.indptr, err_msg=err_msg)
    np.testing.assert_equal(a.indices, b.indices, err_msg=err_msg)    
    np.testing.assert_almost_equal(a.data, b.data, decimal=decimal, err_msg=err_msg)


def test_dissimilarity_matrix(max_grid_size):

    for grid_size in range(max_grid_size + 1):
        # generate an n x n gird of points
        x,y                         =   np.meshgrid( np.arange(grid_size), np.arange(grid_size) )
        cloud                       =   np.column_stack((x.ravel(), y.ravel()))

        # calculate the distance matrix, and format as a csr matrix; we regard these two matrices as "raw data"
        if grid_size > 0:
            dissimilarity_matrix_dense            =   sklearn.metrics.pairwise_distances(cloud)
            dissimilarity_matrix_dense            =   np.maximum( dissimilarity_matrix_dense, dissimilarity_matrix_dense.T) # ensure the matrix is symmetric
        else:
            dissimilarity_matrix_dense            =   np.zeros((0,0))
        dissimilarity_matrix_csr                  =   scipy.sparse.csr_matrix(dissimilarity_matrix_dense)
        dissimilarity_matrix_csr.setdiag(0)

        # format the raw data from the cloud, dense, and csr matrices into "cleaned" csr matrices
        formatted_from_cloud        =   matrix_from_cloud_slow(          cloud       , dissimilarity_max=inf    )    
        formatted_from_cloud_nn     =   matrix_from_cloud(   cloud       , dissimilarity_max=inf    )
        formatted_from_csr          =   matrix_from_csr(            dissimilarity_matrix_csr  , dissimilarity_max=inf    )    
        formatted_from_dense        =   matrix_from_dense(          dissimilarity_matrix_dense, dissimilarity_max=inf    )

        # check that these matrices are all approximately equal
        matrices = [formatted_from_cloud, formatted_from_cloud_nn, formatted_from_csr, formatted_from_dense ]
        for (inda, a), (indb, b) in itertools.product( enumerate(matrices), enumerate(matrices) ):
            assert_almost_equal_csr( a, b, decimal=10, err_msg=f"inda: {inda}, indb: {indb}" )

        # generate a list of distance upper bounds, and check that they are similar to one another       

        radius_cloud_slow           =   enclosing_from_cloud_slow(cloud)
        radius_cloud                =   enclosing_from_cloud(cloud)
        radius_csr                  =   enclosing_from_csr(dissimilarity_matrix_csr)
        radius_dense                =   enclosing_from_dense(dissimilarity_matrix_dense)

        np.testing.assert_almost_equal( radius_cloud_slow, radius_cloud,            decimal=10 )
        np.testing.assert_almost_equal( radius_cloud_slow, radius_csr,              decimal=10 )
        np.testing.assert_almost_equal( radius_cloud_slow, radius_dense,            decimal=10 )    

        # generate some sparsified matrices, and check that they are similar
        thresholds      =\
        [ radius_cloud_slow, radius_cloud, radius_csr, radius_dense    ] +\
        [ 0,  1,  2,  3,  4,  5,  inf               ] +\
        [ 0, -1, -2, -3, -4, -5, -inf               ] +\
        list( np.unique(dissimilarity_matrix_dense)               ) +\
        list( dissimilarity_matrix_csr.data                       ) +\
        list( formatted_from_cloud.data             ) +\
        list( formatted_from_cloud.data             ) +\
        list( formatted_from_cloud_nn.data          ) +\
        list( formatted_from_csr.data               ) +\
        list( formatted_from_dense.data             )   

        for threshold in thresholds:
            formatted_from_cloud        =   matrix_from_cloud_slow(     cloud,                      dissimilarity_max=threshold    )    
            formatted_from_cloud_nn     =   matrix_from_cloud(          cloud,                      dissimilarity_max=threshold + 0.00000000001    )
            formatted_from_csr          =   matrix_from_csr(            dissimilarity_matrix_csr  , dissimilarity_max=threshold    )    
            formatted_from_dense        =   matrix_from_dense(          dissimilarity_matrix_dense, dissimilarity_max=threshold    )        

            matrices = [formatted_from_cloud, formatted_from_cloud_nn, formatted_from_csr, formatted_from_dense ]
            for (inda,a), (indb,b) in itertools.product(enumerate(matrices), enumerate(matrices)):
                assert_almost_equal_csr( a, b, decimal=10, err_msg=f"inda: {inda}, indb: {indb}")

    print("test passed")
