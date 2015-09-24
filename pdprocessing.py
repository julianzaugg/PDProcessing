"""
For processing results from PairDistances
"""


import argparse
import os,sys
import shutil


import annotation
import subgroup
import sequence
import distance_matrix as dm
import tree_constructor as tc

__author__ = 'julianzaugg'


SEQUENCES = None
DISTANCE_MATRIX = None

PD_INPUT_FILENAME = None
OUTPUT_DIR = None

def load_data(input_filename):
    global SEQUENCES
    global DISTANCE_MATRIX
    # Create Annotation object
    my_annotation = annotation.Annotation(input_filename)

    # Grab relevant columns
    # sources = my_annotation.get_column("Seq1")
    targets = my_annotation.get_column("Seq2")
    distances = my_annotation.get_column("Time")

    # s1_content = my_annotation.get_column("S1Sequence")
    s2_content = my_annotation.get_column("S2Sequence")

    # Create sequence objects - ordered
    visited = set()
    seqs = []
    for seq_name, seq_content in zip(targets[::-1], s2_content[::-1]):
        if seq_name not in visited:
            seqs.append(sequence.Sequence(name=seq_name, sequence=seq_content))
            visited.add(seq_name)
        elif seq_name in visited:
            break
    seqs = seqs[::-1]
    SEQUENCES = seqs

    #Create Distance matrix
    DISTANCE_MATRIX = _dm_from_pd_output(seqs, distances)

def create_tree():
    global TREE
    TREE = tc.UPGMA(DISTANCE_MATRIX)
    TREE.calculate_tree()

def _no_dups(alist):
    """
    Return a list with no duplicates with ordered preserved
    """
    out = []
    [out.append(i) for i in alist if not out.count(i)]
    return out

def _dm_from_pd_output(seqs, distances):
    """
    Return a distance matrix object from the sources and distances data from a PairDistances results file
    """
    number_of_entries = len(seqs)
    dist_matrix = dm.DistanceMatrix(seqs)
    cnt = 0
    for i in xrange(number_of_entries):
        for j in range(i):
            dist_matrix[i][j] = dist_matrix[j][i] = distances[cnt]
            cnt += 1
        cnt += 1
    return dist_matrix

def is_valid_file(parser, arg):
    """Check if arg is a valid file that already exists on the file
       system.
    """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg

def _dir_check(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def _make_subfolders(out_dir, clean_dirs=True):
    if clean_dirs:
        shutil.rmtree(out_dir + "/Fastas/")
        shutil.rmtree(out_dir + "/Alignments/")
        shutil.rmtree(out_dir + "/Logos/")
    if not os.path.exists(out_dir + "/Fastas/"):
        os.makedirs(out_dir + "/Fastas/")
    if not os.path.exists(out_dir + "/Alignments/"):
        os.makedirs(out_dir + "/Alignments/")
    if not os.path.exists(out_dir + "/Logos/"):
        os.makedirs(out_dir + "/Logos/")


def _get_group_sequences():
    sg = subgroup.SubGroup(SEQUENCES, OUTPUT_DIR + "tree.txt")
    groups = sg.find_sub_families()
    new_dict = dict()
    for name, names_in_group in groups.items():
        new_dict[name] = [seq for seq in SEQUENCES if seq.name in names_in_group]
    return sg, new_dict


def process_args(input_parser, input_args):
    global PD_INPUT_FILENAME
    global OUTPUT_DIR
    global DISTANCE_MATRIX
    global SEQUENCES
    is_valid_file(input_parser, input_args.input)
    PD_INPUT_FILENAME = input_args.input
    _dir_check(input_args.output)
    OUTPUT_DIR = input_args.output
    _make_subfolders(OUTPUT_DIR)
    # if input_args.distancematrix and input_args.tree:
    #     DISTANCE_MATRIX = dm.read_phylip(input_args.distancematrix)
    #     TREE
    # else:
    print "Loading data and creating distance matrix"
    load_data(PD_INPUT_FILENAME)
    if input_args.limitlength:
        print "Extracting partial distance matrix using length limits"
        _ll = input_args.limitlength.split(',')
        assert len(_ll) == 2, "Incorrectly specified range values"
        min_length, max_length = int(_ll[0]), int(_ll[1])
        exclude_idxs = []
        for i in xrange(len(SEQUENCES)):
            s_length = len(SEQUENCES[i])
            if s_length < min_length or s_length > max_length:
                exclude_idxs.append(i)
        assert len(exclude_idxs) != len(SEQUENCES), "All sequences were excluded, change length limit range"
        DISTANCE_MATRIX = dm.DistanceMatrix.get_partial_dmatrix(DISTANCE_MATRIX, exclude_idxs)
    print "Writing Distance matrix"
    dm.write_phylip(OUTPUT_DIR + "distance_matrix.txt", DISTANCE_MATRIX)
    print "Creating Tree"
    create_tree()
    print "Writing Tree"
    TREE.write_tree(OUTPUT_DIR + "tree.txt")
    print "Identifying subgroups"
    _temp = _get_group_sequences()
    sg_object, groups = _temp[0], _temp[1]
    print "Writing group annotations"
    sg_object.write_annotations(OUTPUT_DIR + "group_annotations.txt")
    print "Writing group fastas"
    sequence.write_fasta_files(OUTPUT_DIR + "Fastas/", groups)
    print "Aligning sequence groups and writing clustal files"
    sequence.align_fasta_and_write(OUTPUT_DIR + "Fastas/", OUTPUT_DIR + "Alignments/")
    print "Generating Logos"
    sequence.generate_logos(OUTPUT_DIR + "Alignments/", OUTPUT_DIR + "Logos/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Description to be added')
    parser.add_argument('-i', '--input', help='Input results file from PairDistances', required=True)
    parser.add_argument('-o', '--output', help='Output location for all results', required=False, default="./")
    parser.add_argument('-ll', '--limitlength', help='Limit processing to sequences of certain length;'
                                                     'Must be in format of "min,max"', required=False)

    args = parser.parse_args()
    process_args(parser, args)
