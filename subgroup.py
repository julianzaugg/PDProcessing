"""
Given a tree this script will identify the super-families using the algorithm in
Krause et. al (2005) - "Large Scale hierarchical clustering of protein sequences"

require - Biopython
"""
__author__ = 'julianzaugg'

# TODO provide option to store all possible groups for each sequence, rather than largest inclusive
# TODO Investigate improvement of algorithm to better identify subgroups. Currently not separating accurately enough.

from Bio import Phylo


class SubGroup:

    def __init__(self, sequences, tree_filename):
        self.seqs = sequences
        self.seqs_dict = dict([(s.name, s) for s in sequences])
        self.tree = self._load_tree(tree_filename)
        self.sub_families = None

    def _load_tree(self, tree_filename):
        return Phylo.read(tree_filename, 'newick')

    def get_parent(self, tree, child_clade):
        node_path = tree.get_path(child_clade)
        return node_path[-2]

    def all_parents(self, tree):
        parents = {}
        for clade in tree.find_clades(order='level'):
            for child in clade:
                parents[child] = clade
        return parents

    # FIXME very ugly method. Efficiency probably could be improved
    def find_sub_families(self, skip_grouped=True):
        """
        Find the super-families in the full tree (or hopefully in our case, subfamilies).
        Based on algorithm presented by Krause et. al (2005)

        If skip_grouped = True, will ignore sequences that are in the super family of another
        sequence.
        """
        already_grouped = set()
        families = {}
        leaves = self.tree.get_terminals()
        parents = self.all_parents(self.tree)
        group = 0
        for leaf in leaves:
            if skip_grouped:
                if leaf.name in already_grouped:
                    continue
            # print "Finding subfamily : " + str(group)
            q = leaf
            I = 0                                               # Largest ratio value stored here
            subfamily = leaf
            while q != self.tree.root:
                parent = parents[q]
                p_count = parent.count_terminals()
                leaf_count = q.count_terminals()
                J = float(p_count - leaf_count) / leaf_count
                if J > I:
                    I = J
                    subfamily = q
                q = parent
            sf_names = [t.name for t in subfamily.get_terminals()]
            families["Group" + str(group)] = (sf_names, subfamily)
            already_grouped = already_grouped.union(set(sf_names))
            group += 1
        # Resolve families by keeping only the largest, i.e., if a subtree is contained within another, remove it
        remove_these = set()
        for family in families:
            f_tree = families[family][1]
            f_tree_size = len(families[family][0])
            for other_family in families:
                of_tree = families[other_family][1]
                of_tree_size = len(families[other_family][0])
                # If the current family subtree includes the `other' family tree, then we don't consider the
                # other family as an individual group and mark it for removal
                if f_tree.is_parent_of(of_tree) and of_tree_size < f_tree_size:
                    remove_these.add(other_family)
        new_families = dict()
        cnt = 0
        for k, v in families.items():
            if k not in remove_these:
                new_name = "Group%i" % cnt
                #Remove [0] to assign a tuple, where [0] is the list of seq names and [1] is the clade object
                new_families[new_name] = v[0]
                cnt += 1
        self.sub_families = new_families
        return new_families

    def change_leaf_names_to_sequence(self, mapping):
        """
        Change the leaf names to the actual sequence
        mapping should be a dictionary with [current_name : sequence] (sequence can also be whatever you want)
        """
        leaves = self.tree.get_terminals()
        for leaf in leaves:
            leaf.name = mapping[leaf.name]

    def write_annotations(self, filename):
        """
        Write an annotation file containing columns for name and the identified subgroup
        :param filename:
        :return:
        """
        assert self.sub_families, "Call 'find_sub_families' first"
        with open(filename, 'w') as fh:
            print >> fh, "Name\tGroup\tSequence"
            for group, results in self.sub_families.items():
                for name in results:
                    print >> fh, "%s\t%s\t%s" % (name, group, self.seqs_dict[name].sequence)

if __name__ == "__main__":
    pass
