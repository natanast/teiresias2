# from memory_profiler import profile
from tqdm import tqdm
import numpy as np
from sklearn import preprocessing

from src.utils.utils import timeit

# fp = open("memory_profiler.log", "w+")

def max_similarity(max_len: int, weight: float, multiplier: float, threshold: int = 105):
    # Maximum similarity for sequence length < 150 is: num_of_aas + ((num_of_aas - 1) * weight)
    _max_similarity_1 = threshold - 1 + (threshold - 2) * weight

    # Maximum similarity for sequence length >= 150 is: num_of_aas * weight + ((num_of_aas - 1) * weight)
    _max_similarity_2 = (max_len - threshold + 1) * multiplier + ((max_len - threshold + 1 ) * weight)

    _sum = _max_similarity_1 + _max_similarity_2

    return _sum


def similarity(seq1: str, seq2: str, weight: float, multiplier: float, penalty: float, threshold: int, gap_threshold: int = 15 ):
    aa_clusters = [
        {"I", "L", "V", "A"},
        {"H", "K", "R"},
        {"M", "C"},
        {"T", "S"},
        {"E", "D"},
        {"Q", "N"},
    ]

    max_concurrencies = 0
    _similarity = 0
    for pos, (aa1, aa2) in enumerate(zip(seq1, seq2), start=1):
        max_concurrencies += 1
        class_multiplier = 1
        
        if pos < gap_threshold:
            if aa1 == ".":
                aa1 = aa2
            
            if aa2 == ".":
                aa2 = aa1
    
        if pos < threshold:
            _current_multiplier = 1
            reward = weight
        else:
            _current_multiplier = multiplier
            reward = weight
            # reward = 1
            
        if aa1 == aa2:   
            if max_concurrencies == 1:
                _similarity += _current_multiplier
            else:
                _similarity += _current_multiplier + reward
        else:
            for cluster in aa_clusters:
                if aa1 in cluster and aa2 in cluster:
                    class_multiplier = 2 / 3
                    break

            if class_multiplier == 2 / 3:
                if max_concurrencies == 1:
                    _similarity += class_multiplier * _current_multiplier
                else:
                    _similarity += (class_multiplier + reward) * _current_multiplier
            else:
                max_concurrencies = 0
                _similarity -= penalty

    _max_similarity = max_similarity(len(seq1), weight, multiplier, threshold)
    # print("Similarity is: ", _similarity)
    # print("Max similarity is: ", _max_similarity)

    _similarity = _similarity / _max_similarity
    return max(min(_similarity,1),0) #np.clip(_similarity, 0, 1)


def distance(seq1: str, seq2: str, weight, multiplier, penalty, threshold, gap_threshold):
    return 1 - similarity(seq1, seq2, weight, multiplier, penalty, threshold, gap_threshold)


def construct_neighbor_matrix(lst):
    length = len(lst) + 1
    distance_matrix = np.zeros((length, length))
    for i, el in enumerate(lst):
        size = len(el)
        for j in range(size):
            distance_matrix[i][j + i + 1] = el[j]
            distance_matrix[i + j + 1][i] = el[j]
    return distance_matrix


@timeit
def distance_matrix(aa_seqs: dict, weight, multiplier, penalty, threshold, gap_threshold):
    total_sequences = len(aa_seqs)
    _distance_matrix = np.zeros((total_sequences, total_sequences))
    #n = int(total_sequences^2 - total_sequences) / 2
    #similarity_matrix = []

    sequences = list(aa_seqs.values())
    snames = list(aa_seqs.keys())

    for idx1, (name1, seq1) in enumerate(
        tqdm(zip(snames, sequences), total=total_sequences, desc="Calculating distance")
    ):
        rest_sequences = sequences[idx1 + 1 :]
        rest_snames = snames[idx1 + 1 :]
        for idx2, (name2, seq2) in enumerate(zip(rest_snames, rest_sequences)):
            distance_score = distance(seq1, seq2, weight, multiplier, penalty, threshold, gap_threshold)
            #similarity_matrix.append(similarity_score)
            _distance_matrix[idx1][idx2 + idx1 + 1] = distance_score
            _distance_matrix[idx1 + idx2 + 1][idx1] = distance_score

    #print(similarity_matrix)
    #normalized_sim_matrix = preprocessing.normalize([similarity_matrix])
    #print(normalized_sim_matrix)
    return _distance_matrix


def max_distance_matrix(distance_matrix: np.ndarray):
    return np.max(distance_matrix)


def min_distance_matrix(distance_matrix: np.ndarray):
    # Remove diagonal
    masked = distance_matrix[~np.eye(distance_matrix.shape[0], dtype=bool)].reshape(
        distance_matrix.shape[0], -1
    )
    return np.min(masked)
