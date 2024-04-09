
import os
import Bio.Align
import numpy as np
import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.phylo as phylo
import biotite.sequence.graphics as graphics
from Bio.Align import substitution_matrices

aligner = Bio.Align.PairwiseAligner(match_score = 1.0)
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
print(aligner.substitution_matrix['Q','Q'])


output = "D:\Documents\MSc\ptyxiaki\\thesis\output.fa"


def read_file(file_path):
    aa_seq = {}
    max_len = 0
    with open(file_path, "r") as f:
        for idx, line in enumerate(f.readlines()):
            if idx % 2 == 0:
                name = line.strip()
            else:
                aa_seq[name] = line.replace('*','').strip()
                if len(aa_seq[name])  > max_len :
                    max_len = len(aa_seq[name])
    return aa_seq, max_len

def seq_padding (aa_seqs, max_len):
    for name, seq in aa_seqs.items():
        aa_seqs[name]=seq.ljust(max_len, ".")
        
    return aa_seqs

def calc_genetic_distance(alignment):

    distance_matrix = []
    for i in range(len(alignment)):

        row = []
        for j in range(i+1):

            seq1 = str(alignment[i])
            seq2 = str(alignment[j])

            # Use Blosum62 matrix to calculate similarity score
            matrix = align.SubstitutionMatrix.std_protein_matrix()
            matrix = align.SubstitutionMatrix(seq.Alphabet(matrix.get_alphabet1().get_symbols()[:-4]), seq.Alphabet(matrix.get_alphabet2().get_symbols()[:-4]),
            matrix.score_matrix()[:-4, :-4])

            #print(matrix)

            #print(type(matrix))

            gap_open = -10
            gap_extend = -0.5
            #print("seq2: ",seq2)
            #print("seq1: ",seq1)
            
            #score = aligner.score(seq1, seq2)
            #print(f"The score of the  alignment is {score}")
            #alignments = PairwiseAligner.align(seq1, seq2, matrix, gap_open, gap_extend)
            alignments = aligner.align(seq1, seq2)

            # Calculate genetic distance as 1 - similarity score
            distance = 1 - (alignments[0].score / len(seq1))

            row.append(distance)

        distance_matrix.append(row)

    return distance_matrix



if __name__ == "__main__":

    aa_seqs, max_len = read_file(output)
    #padded_seq = seq_padding(aa_seqs, max_len)
    my_array = np.array(list(aa_seqs.values()))
    print(len(my_array))
    calc_genetic_distance(my_array)
    #print(calc_genetic_distance(my_array))

    #print(padded_seq)
    #print(calc_genetic_distance())