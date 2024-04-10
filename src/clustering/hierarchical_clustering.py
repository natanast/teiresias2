import os
from datetime import datetime
import multiprocessing
import biotite.sequence.graphics as graphics
import biotite.sequence.phylo.upgma as upgma
import biotite.sequence as seq
import matplotlib.pyplot as plt
# import array
import numpy as np
from time import time

# from src.utils.utils import timeit

# @timeit
def get_tree(dist_matrix: np.ndarray, save: bool = True):
    start = time()

    tree = upgma(dist_matrix)
    if save:
        current_datetime = datetime.now()
        datetime_str = current_datetime.strftime("%Y%m%d%H%M%S")
        file_path = os.path.join(
            os.getcwd(), "results", "trees", "tree_" + datetime_str + ".nwk"
        )
        save_tree_structure(tree, file_path)
    end = time()
    print(f"get_tree took {end - start} seconds")

    return tree, file_path


def save_tree_structure(tree, file_path):
    with open(file_path, "w") as f:
        f.write(tree.to_newick(include_distance=True))


def visualize_tree(tree, aa_seqs: dict):
    fig = plt.figure(figsize=(10.0, 5.0))
    ax = fig.add_subplot(111)
    graphics.plot_dendrogram(ax, tree, orientation="left", labels=list(aa_seqs.keys()))
    # plt.ylabel("Sequences",  labelpad=1, font_size = 4)

    # ax.set_xlabel("Sequences",labelpad=5)
    plt.xticks(rotation=90)
    # ax.set_xticklabels(ax.get_xticks(), rotation = 90)

    # ax.yaxis.grid(color="lightgray")
    plt.savefig("mygraph.png", bbox_inches="tight", dpi=300, pad_inches=0.5)
