from typing import List, Optional
from pandas import DataFrame



class FactoredBoundaryMatrixDowker:
    """
    The factored boundary matrix of a Dowker complex
    :param dowker_simplices a list of softed-in-ascending-order lists of integers
    :param max_homology_dimension the maximum dimension in which we want to compute homology
    """

    def __init__(self, dowker_simplices: List[List[int]], max_homology_dimension: int): 
        """
        Initializes a new FactoredBoundaryMatrixDowker instance with a list of (sorted) lists of integers.

        Args:
            dowker_simplices: A list of (sorted in strictly ascending order) lists of integers to use for initialization.
        """

    def homology(self) -> DataFrame:
        """
        Returns a Pandas DataFrame with information about homology, betti numbers, and cycle representatives
        """

    # def jordan_column_for_simplex( keymaj: List ) -> DataFrame:
    #     """
    #     Obtain a column of the Jordan basis associated with the U-match factorization

    #     :param keymaj: a list of integers
        
    #     :return L: a list of tuples `( s, a, b )`, where `s` is a simplex, and `a/b` is the coefficient of the simplex
    #     """


class FactoredBoundaryMatrixVr:
    """
    The factored boundary matrix of a filtered Vietoris Rips complex
    
    This object is uniquely determined by three user-defined parameters:
    - a dissimilarity_matrix
    - a maximum dissimilarity threshold
    - the maximum dimension in which we want to compute homology
    """

    def __init__(self,  dissimilarity_matrix: List[List[float]], dissimilarity_max: Optional[float],  homology_dimension_max:  int, ):             
        """
        Initializes a new FactoredBoundaryMatrixVr instance

        Args:
            dissimilarity_matrix: a symmetric matrix (such as a distance matrix) represented in list-of-list format
            dissimilarity_max: we only consttruct simplices with diameter diameter or smaller
            homology_dimension_max: the maximum dimension to compute homology
        """

    def homology(self) -> DataFrame:
        """
        Returns a Pandas DataFrame with information about homology, betti numbers, and cycle representatives
        """
 

    # def jordan_column_for_simplex( keymaj: List ) -> DataFrame:
    #     """
    #     Obtain a column of the Jordan basis associated with the U-match factorization

    #     :param keymaj: a list of integers
        
    #     :return L: a list of tuples `( s, a, b )`, where `s` is a simplex, and `a/b` is the coefficient of the simplex
    #     """        