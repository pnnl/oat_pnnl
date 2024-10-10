
from cmath import inf



def max_finite_value( values ):
    """
    Find the maximum finite value less than +infinity, in a list. Returns None if the list does not contain a
    a value < infinity.

    :param values: an iterable of floats
    """
    finite_points = [ v for v in values if v < inf ]
    if len(finite_points)==0:
        return None
    else:
        return max(  finite_points )




# def max_finite_endpoint( birth_death_pairs ):
#     """
#     Find the maximum finite endpoint of any bar in the barcode.
#     If the barcode is empty then return None.

#     :param birth_death_pairs: an iterable of (birth,death) pairs
#     """
#     finite_points = [ x[p] for x in birth_death_pairs for p in [0,1] if x[p] < np.inf ]
#     if len(finite_points)==0:
#         return None
#     else:
#         return max(  finite_points )
