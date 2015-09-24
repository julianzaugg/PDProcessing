

__author__ = 'julianzaugg'

import numpy as np
from itertools import chain

class DistanceMatrix:
    """
    Stores a distance matrix. Allows printing of matrix in PHYLIP formatted distance matrix.

    """

    def __init__(self, seqs, force_use_name=False):
        """
        :param seqs: Sequences associated with this distance matrix. Assumes ordered ot match row/column indexing.
        :param force_use_name: force the use of the sequence names (truncated)
        """
        self._seqs = seqs
        self.force_name = force_use_name
        self.n_seqs = len(seqs)
        self.dist_matrix = np.zeros((self.n_seqs, self.n_seqs))
        self.names = [s.name for s in seqs]
        self.name_map = dict()

    def __str__(self):
        """
        Return a string representation of the matrix in phylip format.
        Each row is assigned an arbitrary name, e.g, Seq0, Seq1..etc, but the actual names can be accessed in
        self.name_map. This is done in case supplied names create duplicates when truncated at 10 characters (which
        is required for PHYLIP formatted matrices).
        """
        out_string = str(self.n_seqs) + "\n"
        for i in xrange(self.n_seqs):
            if not self.force_name:
                s1_name = "Seq%i" % i
            else:
                s1_name = self._seqs[i].name[:10]
            while len(s1_name) < 10:
                s1_name += " "
            self.name_map[s1_name] = self._seqs[i].name
            row_data_string = "\t".join(map(str, self.dist_matrix[i][:i]))
            row_string = s1_name + "\t" + row_data_string + "\n"
            out_string += row_string
        return out_string

    def __len__(self):
        return len(self.dist_matrix)

    def __setitem__(self, index, new_value):
        self.dist_matrix[index] = new_value

    def __getitem__(self, index):
        """
        Get row from matrix
        :param index: row index
        :return: Return list of distances in respect to row index
        """
        return self.dist_matrix[index]

    def get_name_mapping(self):
        """
        :return: A dictionary containing the mapping from assigned sequence name to the original name
        """
        return self.name_map

    def save_matrix(self, filename):
        """
        Save the PHYLIP matrix to file
        :param filename: Save location filename
        """
        with open(filename, 'w') as fh:
            print >> fh, str(self)

    def get_lower(self):
        lower = []
        for i in range(len(self)):
            row = []
            for j in range(i + 1):
                row.append(self[i][j])
            lower.append(row)
        return lower

    @staticmethod
    def get_partial_dmatrix(current_dmatrix, exclude_idx):
        """
        Return a sub matrix from the full matrix, excluding specified row and column
        :param exclude_idx: row/column index to exclude
        :return: partial matrix
        """
        N = len(current_dmatrix)
        # new_seqs = [current_dmatrix._seqs[i] for i in xrange(N) if i not in exclude_idx]
        new_seqs = []
        rows_idxs = []
        for i in xrange(N):
            if i not in exclude_idx:
                new_seqs.append(current_dmatrix._seqs[i])
                rows_idxs.append([i])
        new_partial_matrix = DistanceMatrix(new_seqs, force_use_name=current_dmatrix.force_name)
        new_partial_matrix.dist_matrix = current_dmatrix.dist_matrix[rows_idxs, list(chain(*rows_idxs))]
        return new_partial_matrix

def read_phylip(filename):
    """
    Load a phylip formatted distance matrix
    :param filename:
    :return:
    """
    with open(filename, 'r') as fh:
        data = [line.strip().split("\t") for line in fh.readlines()]
        n_seqs = int(data[0])
        new_matrix = DistanceMatrix(n_seqs)
        for l in range(n_seqs):
            new_matrix.names.append(data[l + 1][0])
            for j in range(l):
                new_matrix[l][j] = data[l + 1][j + 1]
    return new_matrix

def write_phylip(filename, distance_matrix):
    """
    Write phylip formatted distance matrix to file
    :param filename:
    :return:
    """
    assert isinstance(distance_matrix, DistanceMatrix), "Provided distance matrix must be of type DistanceMatrix"
    with open(filename, 'w') as fh:
        print >> fh, distance_matrix





if __name__ == "__main__":

    import sequence
    s1 = sequence.Sequence(name = "A", sequence = "A")
    s2 = sequence.Sequence(name = "B", sequence = "DDG")
    s3 = sequence.Sequence(name = "C", sequence = "CEA")
    s4 = sequence.Sequence(name = "C", sequence = "CEA")
    s5 = sequence.Sequence(name = "C", sequence = "CEA")
    dm = DistanceMatrix([s1,s2,s3, s4, s5])
    dm[2] = 1
    dm[4][2] = 3
    print dm.dist_matrix
    out = DistanceMatrix.get_partial_dmatrix(dm, [0,3])
    print out