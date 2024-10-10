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

    def optimize_cycle(self, birth_simplex: List[int] ) -> dict:
        """        
        Optimize a cycle representative

        Specifically, we employ the "edge loss" method to find a solution `x'` to the problem 

        `minimize Cost(Ax + z)`

        where 
            - `x` is unconstrained
        - `z` is a cycle representative for a (persistent) homology class associated to `birth_simplex`
        - `A` is a matrix composed of a subset of columns of the Jordna basis
        - `Cost(z)` is the sum of the absolute values of the products `z_s * diameter(s)`.

        # Arguments

        - The `birth_simplex` of a cycle represenative `z` for a bar `b` in persistent homology.
        - The `constraint` type for the problem. The optimization procedure works by adding linear
        combinations of column vectors from the Jordan basis matrix computed in the factorization.
        This argument controls which columns are available for the combination.
          - (default) **"preserve PH basis"** adds cycles which appear strictly before `birth_simplex`
            in the lexicographic ordering on filtered simplex (by filtration, then breaking ties by
            lexicographic order on simplices) and die no later than `birth_simplex`.  **Note** this is
            almost the same as the problem described in [Escolar and Hiraoka, Optimal Cycles for Persistent Homology Via Linear Programming](https://link.springer.com/chapter/10.1007/978-4-431-55420-2_5)
            except that we can include essential cycles, if `birth_simplex` represents an essential class. 
          - **"preserve PH basis (once)"** adds cycles which (i) are distince from the one we want to optimize, and
            (ii) appear (respectively, disappear) no later than the cycle of `birth_simplex`.  This is a looser
            requirement than "preserve PH basis", and may therefore produce a tighter cycle.  Note,
            however, that if we perform this optimization on two or more persistent homology classes in a
            basis of cycle representatives for persistent homology, then the result may not be a
            persistent homology basis.
          - **"preserve homology class"** adds every boundary vector
          - "preserve homology calss (once)" adds every cycle except the one represented by `birth_simplex`



        # Returns

        - The vectors `b`, `x`, and `y`
          - We separate `x` into two components: one made up of codimension-1 simplices (labeled "difference in bounding chains"), and one made up of codimension-0 simplices (labeled "difference in essential cycles")
        - The objective values of the initial and optimized cycles
        - The number of nonzero entries in the initial and optimized cycles
    
    # Related
    
    See
    
    - [Escolar and Hiraoka, Optimal Cycles for Persistent Homology Via Linear Programming](https://link.springer.com/chapter/10.1007/978-4-431-55420-2_5)
    - [Obayashi, Tightest representative cycle of a generator in persistent homology](https://epubs.siam.org/doi/10.1137/17M1159439)
    - [Minimal Cycle Representatives in Persistent Homology Using Linear Programming: An Empirical Study With User’s Guide](https://www.frontiersin.org/articles/10.3389/frai.2021.681117/full)
              # def jordan_column_for_simplex( keymaj: List ) -> DataFrame:
    #     
    """