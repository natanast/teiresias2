#python main.py C:\Users\Nina\ptyxiaki\thesis\data\output.fa 1.5
#python main.py C:\Users\Nina\ptyxiaki\thesis\data\3000_seqs.fa 2.0
#python main.py C:\Users\Nina\ptyxiaki\thesis\data\tiny_test.fa 2.0
#python main.py C:\Users\Nina\ptyxiaki\thesis\data\62K.fa 2.0

import os
from datetime import datetime
import argparse
from src.clustering.hierarchical_clustering import get_tree, visualize_tree
from src.metrics.metrics import distance_matrix
from src.utils.cli import parse_argument
from src.utils.utils import read_file, seq_padding, timeit
from src.clustering.tree_2_clusters import divide_tree_into_clusters
import numpy as np
from Bio import Phylo
import phylotreelib as pt
#from src.clustering.hierarchical_clustering import readDistMatrix, upgma, printCluster

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="Input file containing AA sequences.")
    parser.add_argument("multiplier", type=float, default=2, help="Multiplier for similarity calculation.")
    args = parser.parse_args()
    #args = parse_argument()
    if not args.input:
        raise Exception("Please specify an input file")
    aa_seqs, max_len = read_file(args.input)
    padded_seq = seq_padding(aa_seqs, max_len)

    print("Total number of sequences:", len(aa_seqs))
    dist_matrix = distance_matrix(padded_seq, max_len, args.multiplier)
    np.savetxt('UPGMA_Input.txt',dist_matrix,fmt='%.4f')
    print("Distance matrix saved")

    print("Creating tree")
    tree = get_tree(dist_matrix) 
    print("Tree created")
    ##########visualize_tree(tree, aa_seqs)

    # print("Loading UPGMA")
    # dist_matrix = np.loadtxt('UPGMA_Input.txt')
    # print("UPGMA loaded")

    # tree_path = "C:/Users/Nina/ptyxiaki/thesis/results/trees/tree_20231203203735.nwk"
    
    # divide_tree_into_clusters(tree_path, padded_seq, dist_matrix)

    ###########tree = Phylo.read(tree_path, "newick")
    ###########Phylo.draw(tree)
