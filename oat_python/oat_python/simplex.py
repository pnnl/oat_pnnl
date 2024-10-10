
import itertools
import numpy as np






def dnfaces( simplices=[], facedim=0 ):
    """
    Return all dimension-`facedim` faces of the simplices provided

    :param simplices: an iterable of iterables
    :param facedim: dimension of desired faces

    :return faces: a set of tuples representing the collection of dimension-`facedim` faces
    """
    faces = set()
    for simplex in simplices:
        for face in itertools.combinations( simplex, facedim+1 ):
            faces.add( face )
    return faces


def dmsimplicesasrows__dnfacesasrows( dksimplices_as_rows: np.ndarray, facedim = 0, removeduplicatefaces=True ):
    """
    :param dksimplices_as_rows: a numpy array where each row represents a simplex 

    :return a numpy array whose rows represent the dimension-`facedim` faces of the simplices

    NB: the returned array will have duplicate rows if two or more simplices share a face, unless removeduplicatefaces=True
    """
    if not isinstance(dksimplices_as_rows, np.ndarray):
        raise TypeError("Input must be a NumPy array.")
    if dksimplices_as_rows.ndim != 2 or dksimplices_as_rows.shape[1] != 3:
        raise ValueError("Array must have shape (n, 3).")

    if len(dksimplices_as_rows) == 0:
        return []
    
    dksimplices_as_rows     =   np.array(dksimplices_as_rows)
    nvertices               =   len(dksimplices_as_rows[0])
    subset_indices          =   itertools.combinations( range(nvertices), facedim+1 )
    dnfacesasrows           =   [dksimplices_as_rows[:,I] for I in subset_indices]
    dnfacesasrows           =   np.concatenate(dnfacesasrows, axis=0)

    if removeduplicatefaces:
        dnfacesasrows       =   np.unique(dnfacesasrows, axis=1)

    return dnfacesasrows