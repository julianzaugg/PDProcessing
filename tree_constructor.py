__author__ = 'julianzaugg'


from Bio.Phylo import TreeConstruction
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
class UPGMA:

    def __init__(self, d_matrix):
        self.matrix = TreeConstruction._DistanceMatrix(names=d_matrix.names, matrix=d_matrix.get_lower())
        self.tree = None

    def __str__(self):
        return self.tree.format("newick")

    def calculate_tree(self):
        constructor = DistanceTreeConstructor()
        self.tree = constructor.upgma(self.matrix)
        return self.tree

    def write_tree(self, filename):
        with open(filename, "w") as fh:
            print >> fh, self