


import unittest
from . import dissimilarity


class TestCase(unittest.TestCase):
    def test_dissimilarity(self):
        dissimilarity.test_dissimilarity_matrix(max_grid_size=3)




if __name__ == '__main__':
    unittest.main()        


