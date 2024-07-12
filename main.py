import argparse
from src.clustering.hierarchical_clustering import get_tree
from src.metrics.metrics import distance_matrix
from src.utils.utils import read_file, seq_padding
from src.clustering.tree_2_clusters import divide_tree_into_clusters
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, help = "Input file containing AA sequences.")
    parser.add_argument("-m", "--multiplier", type = float, default = 2.0, help = "Multiplier for similarity calculation after threshold aminoacid position")
    parser.add_argument("-w", "--weight", type = float, default= 0.1, help = "Added weight for concurrent matches")
    parser.add_argument("-p", "--penalty", type = float, default = 0.0, help = "Penalty for missmatch")
    parser.add_argument("-t", "--threshold", type = int, default = 105, help = "Amino acid position number after which the match will be multiplied with 'multiplier' argument")
    parser.add_argument("-c", "--cuts", type = int, choices = [4, 8, 12], default = 4, help = "Number of cuts in the hierarchicall tree (4, 8, or 12)")
    parser.add_argument("-o", "--output", type = str, help = "Excel file containing the initial data, where calculated clusters will be written as new columns")
    parser.add_argument("-g", "--gap_threshold", type = int, default = 15, help = "Threshold up to which gaps are considered matches")

    args = parser.parse_args()

    # Extracting arguments from the user input
    input_file = args.input
    w = args.weight
    m = args.multiplier
    p = args.penalty
    t = args.threshold
    c = args.cuts
    o = args.output
    g = args.gap_threshold
    
    if not input_file:
        raise Exception("Please specify an input file")
 
    aa_seqs, max_len = read_file(input_file)
    padded_seq = seq_padding(aa_seqs, max_len)

    print("Total number of sequences:", len(aa_seqs))

    dist_matrix = distance_matrix(padded_seq, w, m, p, t, g)
    np.savetxt('UPGMA_Input.txt',dist_matrix,fmt='%.4f')
    print("Distance matrix saved\n")

    print("Creating tree...")
    tree, tree_path = get_tree(dist_matrix) 
    print("Tree created\n")
    
    divide_tree_into_clusters(tree_path, padded_seq, dist_matrix, c, o)