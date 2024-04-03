import biotite.sequence.graphics as graphics
import biotite.sequence.phylo.upgma as upgma
import matplotlib.pyplot as plt
import numpy as np
# from Bio.Align import substitution_matrices
# from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

# aligner = Bio.Align.PairwiseAligner(match_score=1.0)
# aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
# print(aligner.substitution_matrix["Q", "Q"])

output = "D:\Documents\MSc\ptyxiaki\\thesis\\test.fa"

def read_file(file_path):
    aa_seq = {}
    max_len = 0
    with open(file_path, "r") as f:
        for idx, line in enumerate(f.readlines()):
            if idx % 2 == 0:
                name = line.strip()
            else:
                aa_seq[name] = line.replace("*", "").strip()
                if len(aa_seq[name]) > max_len:
                    max_len = len(aa_seq[name])
    return aa_seq, max_len


def seq_padding(aa_seqs, max_len):
    for name, seq in aa_seqs.items():
        aa_seqs[name] = seq.ljust(max_len, ".")

    return aa_seqs


def construct_neighbor_matrix(lst):
    length = len(lst) + 1
    distance_matrix = np.zeros((length, length))
    for i, el in enumerate(lst):
        size = len(el)
        for j in range(size):
            distance_matrix[i][j + i + 1] = el[j]
            distance_matrix[i + j + 1][i] = el[j]
    return distance_matrix


def compare_aas(aa_seqs, max_len):
    tr_distance_matrix = []

    for i in range(len(aa_seqs) - 1):
        row = []

        aa_dict_values = list(aa_seqs.values())

        for j in range(i + 1, len(aa_seqs)):
            seq1 = str(aa_dict_values[i])
            seq2 = str(aa_dict_values[j])
            # print(seq1)
            # print(seq2)

            k = 1
            similarity = 0

            for aai, aaj in zip(seq1, seq2):
                
            # THELW NA ALLAKSW TO ELSE ME TIMH APO PINAKA BLOSUM? # SIMILARITY NOT DIST
            # distance = sum(
            #     (k+1) if aai == aaj else 0 for k, (aai, aaj) in enumerate(zip(seq1, seq2))
            # ) / float(max_len)
            distance = row.append(distance)
            print(i, j, distance)
        tr_distance_matrix.append(row)

    # print(tr_distance_matrix[1][0])
    distance_matrix = construct_neighbor_matrix(tr_distance_matrix)

    print(distance_matrix)
    return distance_matrix

if __name__ == "__main__":
    aa_seqs, max_len = read_file(output)
    padded_seq = seq_padding(aa_seqs, max_len)
    print(len(aa_seqs))
    # for values in padded_seq.values():
    #    for i in range(len(values)):
    #         print(values[i])

    distMatrix = compare_aas(padded_seq, max_len)
    out_tree = upgma(distMatrix)
    
    #print(out_tree.to_newick(include_distance=True))
    file_path = "tree_output.nwk"

    with open(file_path, "w") as file:
    
        file.write(out_tree.to_newick(include_distance=True))

    fig = plt.figure(figsize=(10.0, 5.0))
    ax = fig.add_subplot(111)

    graphics.plot_dendrogram(
        ax, out_tree, orientation="top", labels=list(aa_seqs.keys())
    )
    #plt.ylabel("Sequences",  labelpad=1, font_size = 4)

    #ax.set_xlabel("Sequences",labelpad=5)
    plt.xticks(rotation='vertical')
    

    #ax.yaxis.grid(color="lightgray")
    plt.show()

