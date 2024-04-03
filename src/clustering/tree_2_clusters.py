import phylotreelib as pt
from phylotreelib import Tree
import numpy as np
import pandas as pd 
import re
import ete3

def extract_similarity_matrix(tree):
    similarities = []
    #print(tree)
    with open('test.txt', 'w') as f :
        f.write(tree.__str__())

    with open('test.txt', 'r') as f :
        lines = f.readlines()

    for line in lines[3:]:
        _line = line.strip()
        elements = _line.split('|')

        if len(elements) <4 :
            break

        idx1 = float(_line.split('|')[1])
        idx2 = float(_line.split('|')[2])
        idx3 = float(_line.split('|')[3])
        similarities.append((idx1, idx2, idx3))

    return np.array(similarities)

def compute_quartiles(similarities_array):
    sorted_similarities = np.sort(similarities_array)
    q1 = np.percentile(sorted_similarities, 25)
    q2 = np.percentile(sorted_similarities, 50)
    q3 = np.percentile(sorted_similarities, 75)

    return q1, q2, q3

def compute_octamores(similarities_array):
    sorted_similarities = np.sort(similarities_array)
    q1 = np.percentile(sorted_similarities, 12.5)
    q2 = np.percentile(sorted_similarities, 25)
    q3 = np.percentile(sorted_similarities, 37.5)
    q4 = np.percentile(sorted_similarities, 50)
    q5 = np.percentile(sorted_similarities, 62.5)
    q6 = np.percentile(sorted_similarities, 75)
    q7 = np.percentile(sorted_similarities, 87.5)

    return q1, q2, q3, q4, q5, q6, q7

def compute_dozens(similarities_array):
    sorted_similarities = np.sort(similarities_array)
    q1 = np.percentile(sorted_similarities, 8.33)
    q2 = np.percentile(sorted_similarities, 16.66)
    q3 = np.percentile(sorted_similarities, 25)
    q4 = np.percentile(sorted_similarities, 33.33)
    q5 = np.percentile(sorted_similarities, 41.66)
    q6 = np.percentile(sorted_similarities, 50)
    q7 = np.percentile(sorted_similarities, 58.33)
    q8 = np.percentile(sorted_similarities, 66.66)
    q9 = np.percentile(sorted_similarities, 75)
    q10 = np.percentile(sorted_similarities, 83.33)
    q11 = np.percentile(sorted_similarities, 91.66)

    return q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11

def get_value_from_match(array, idx1, idx2):
    distance = array.astype(float)[:, 2]
    target= np.array([idx1,idx2]).astype(int)
    search_array = array[:,  :-1].astype(int)
    #print(target)
    for i, row in enumerate(search_array):
        #print("DIST OUTSIDE", distance[i])
        if np.array_equal(row, target) or np.array_equal(row, target[::-1]):
            #print(distance[i])
            return distance[i]
    return -1

def compute_cluster_similarity(cluster, distance_matrix):
    num_seqs = len(cluster)
    #print("NUM_SEQS:", num_seqs)
    total_similarity = 0.0
    sim = []

    if num_seqs >= 2 :
        for i in range(num_seqs):
            for j in range(i+1, num_seqs):
                seq_i = cluster[i]
                seq_j = cluster[j]
                # print("########################################################################")
                # print(seq_i)
                #PREVIOUS VERSION: 
                similarity = 1 - distance_matrix[int(seq_i)][int(seq_j)]
                #similarity =  similarity_matrix[int(seq_i)][int(seq_j)]
                #print("SIMILARITY: ",   similarity)
                sim.append(similarity)
        
        total_similarity = sum(sim) / len(sim)
    #print("TOTAL SIMILARITY: ",   total_similarity)

    return total_similarity

def compute_outer_cluster_similarity(cluster, distance_matrix, all_elements):
    num_seqs = len(cluster)
    total_outer_similarity = 0.0
    outer_sim = []
    # all_elements = [i for i in range(len(distance_matrix))]

    if num_seqs >= 1:
        for i in range(num_seqs):
            seq_i = cluster[i]

            for j in range(len(all_elements)):
                if j not in cluster:  # Compare with elements outside the cluster
                    seq_j = all_elements[j]
                    similarity = 1 - distance_matrix[int(seq_i)][int(seq_j)]
                    outer_sim.append(similarity)

        total_outer_similarity = sum(outer_sim) / len(outer_sim) if len(outer_sim) > 0 else 0.0

    return total_outer_similarity

def get_tree_length (tree_path):
    tree = ete3.Tree(tree_path)
    farthest, dist = tree.get_farthest_leaf()

    return dist

def divide_tree_into_clusters (tree_path, fasta_seqs, dist_matrix): 

    with pt.Newicktreefile(tree_path) as treefile:
        tree = treefile.readtree()

    with open(tree_path, 'r') as file:
        tree_str = file.read()
    branch_len_values = [float(match.group(1)) for match in re.finditer(r':([-+]?\d*\.\d+([eE][-+]?\d+)?)', tree_str)]

    tree_len = get_tree_length(tree_path)

    print("\nTree length is: \n", tree_len)

    unique_float_valuesss = list(set(branch_len_values))
    
    # print(max(unique_float_valuesss))
    unique_float_values = [tree_len - value for value in unique_float_valuesss]
    # print(unique_float_values)

    seqs_similarity_nums_col = unique_float_values
    seqs_similarity_matrix = extract_similarity_matrix(tree)

    #seqs_similarity_nums_col = seqs_similarity_matrix[:, 2]
    
    # q1, q2, q3 = compute_quartiles(seqs_similarity_nums_col)
    q1, q2, q3, q4, q5, q6, q7 = compute_octamores(seqs_similarity_nums_col)
    # q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11 = compute_dozens(seqs_similarity_nums_col)

    # h1 = float ((min(seqs_similarity_nums_col) + q1) / 2)
    # h2 = float ( (q1 + q2) / 2 )
    # h3 = float ( (q2 + q3) / 2 )
    # h4 = float ( (q3 + max(seqs_similarity_nums_col)) / 2)

    h1 = float ((min(seqs_similarity_nums_col) + q1) / 2)
    h2 = float ( (q1 + q2) / 2 )
    h3 = float ( (q2 + q3) / 2 )
    h4 = float ( (q3 + q4) / 2)
    h5 = float ( (q4 + q5) / 2)
    h6 = float ( (q5 + q6) / 2 )
    h7 = float ( (q6 + q7) / 2 )
    h8 = float ( (q7 + max(seqs_similarity_nums_col)) / 2)

    # h1 = float ((min(seqs_similarity_nums_col) + q1) / 2)
    # h2 = float ( (q1 + q2) / 2 )
    # h3 = float ( (q2 + q3) / 2 )
    # h4 = float ( (q3 + q4) / 2)
    # h5 = float ( (q4 + q5) / 2)
    # h6 = float ( (q5 + q6) / 2 )
    # h7 = float ( (q6 + q7) / 2 )
    # h8 = float ( (q7 + q8) / 2)
    # h9 = float ( (q8 + q9) / 2)
    # h10 = float ( (q9 + q10) / 2 )
    # h11 = float ( (q10 + q11) / 2 )
    # h12 = float ( (q11 + max(seqs_similarity_nums_col)) / 2)

    # h1 = float(tree_len/8)
    # h2 = float( h1 + tree_len/4 )
    # h3 = float( h2 + tree_len/4 )
    # h4 = float( h3 + tree_len/4 )

    # h1 = float (tree_len/16)
    # h2 = float ( h1 + tree_len/8 )
    # h3 = float ( h2 + tree_len/8 )
    # h4 = float ( h3 + tree_len/8 )
    # h5 = float ( h4 + tree_len/8 )
    # h6 = float ( h5 + tree_len/8 )
    # h7 = float ( h6 + tree_len/8 )
    # h8 = float ( h7 + tree_len/8 )

    # h1 = float (tree_len/24)
    # h2 = float ( h1 + tree_len/12 )
    # h3 = float ( h2 + tree_len/12 )
    # h4 = float ( h3 + tree_len/12 )
    # h5 = float ( h4 + tree_len/12 )
    # h6 = float ( h5 + tree_len/12 )
    # h7 = float ( h6 + tree_len/12 )
    # h8 = float ( h7 + tree_len/12 )
    # h9 = float ( h8 + tree_len/12 )
    # h10 = float ( h9 + tree_len/12 )
    # h11 = float ( h10 + tree_len/12 )
    # h12 = float ( h11 + tree_len/12 )

    # for height in [h1, h2, h3, h4]:
    #     print(height)

    # for height in [h1, h2, h3, h4, h5, h6, h7, h8]:
    #     print(height)

    # heights = {"h1": h1,"h2": h2,"h3": h3, "h4": h4}
    heights = {"h1": h1,"h2": h2,"h3": h3, "h4": h4, "h5": h5,"h6": h6,"h7": h7, "h8": h8}
    # heights = {"h1": h1,"h2": h2,"h3": h3, "h4": h4, "h5": h5,"h6": h6,"h7": h7, "h8": h8, "h9": h9,"h10": h10,"h11": h11, "h12": h12}

    # for height in ['h1', 'h2', 'h3', 'h4']:    
    for height in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'h7', 'h8']:
    # for height in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'h7', 'h8', 'h9', 'h10', 'h11', 'h12']:
        if height == 'h1':
            cutted_tree = Tree.cluster_cut(tree, h1)
        elif height == 'h2':
            cutted_tree = Tree.cluster_cut(tree, h2)
        elif height == 'h3':
            cutted_tree = Tree.cluster_cut(tree, h3)
        elif height == 'h4':
            cutted_tree = Tree.cluster_cut(tree, h4)
        elif height == 'h5':
            cutted_tree = Tree.cluster_cut(tree, h5)
        elif height == 'h6':
            cutted_tree = Tree.cluster_cut(tree, h6)
        elif height == 'h7':
            cutted_tree = Tree.cluster_cut(tree, h7)
        elif height == 'h8':
            cutted_tree = Tree.cluster_cut(tree, h8)
        # elif height == 'h9':
        #     cutted_tree = Tree.cluster_cut(tree, h9)
        # elif height == 'h10':
        #     cutted_tree = Tree.cluster_cut(tree, h10)
        # elif height == 'h11':
        #     cutted_tree = Tree.cluster_cut(tree, h11)
        # elif height == 'h12':
        #     cutted_tree = Tree.cluster_cut(tree, h12)

        clusters = cutted_tree[0]  # No. of clusters
        # print(clusters)
        
        single_element_clusters = [cluster for cluster in clusters if len(cluster) == 1]
        multiple_element_clusters = [cluster for cluster in clusters if len(cluster) > 1]

        # print(single_element_clusters)
        
        Dict = {}
        excel_file = "C:/Users/Nina/ptyxiaki/thesis/data/CLL-DB-data-aligned_random_3000.xlsx"
        with open(f"C:/Users/Nina/ptyxiaki/thesis/results/clusters/final_clusters_{height}.fa", 'w') as f:
            
            for cluster_index, cluster in enumerate(multiple_element_clusters):
                cluster_similarity  = compute_cluster_similarity(list(cluster), dist_matrix)
                all_elements = [i for i in range(len(dist_matrix))]
                dissimilarity = compute_outer_cluster_similarity(list(cluster), dist_matrix, all_elements)
                
                f.write("Cluster {}, Num of seqs: {}, Similarity degree: {} Dissimilarity: {}\n".format(cluster_index + 1, len(cluster),  cluster_similarity, dissimilarity))
                 
                for entry in cluster:
                    key = list(fasta_seqs.keys())[int(entry)]
                    Dict[key] = cluster_index + 1
                    #value = fasta_seqs[key]
                    f.write("{}\n".format(key))
                    #f.write("{}\n".format(value))
                f.write("\n")

            for cluster_index, cluster in enumerate(single_element_clusters):
                cluster_similarity  = "-"
                
                f.write("Cluster {}, Num of seqs: {}, Similarity degree: {} \n".format(cluster_index + len(multiple_element_clusters) + 1, len(cluster),  cluster_similarity))
                    
                for entry in cluster:
                    key = list(fasta_seqs.keys())[int(entry)]
                    Dict[key] = 0
                    #value = fasta_seqs[key]
                    f.write("{}\n".format(key))
                    #f.write("{}\n".format(value))
                f.write("\n")
                
        # sorted_dict = dict(sorted(Dict.items()))
        # ##print(sorted_dict)

        # values_list = [value for key, value in sorted_dict.items()]
        # ##print(values_list)

        # df = pd.read_excel(excel_file)
        # column_name = "Predicted Cluster for h = {}".format(heights[height])
        # df[column_name] = values_list
        # df.to_excel(excel_file, index=False)

    # for height in ['h1', 'h2', 'h3', 'h4']:
    for height in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'h7', 'h8']:
    # for height in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'h7', 'h8', 'h9', 'h10', 'h11', 'h12']:
        with open(f"C:/Users/Nina/ptyxiaki/thesis/results/clusters/final_clusters_{height}.fa", "r") as file:
            clusters_content = file.readlines()

        excel_data = pd.read_excel(excel_file)

        cluster_mapping = {}
        current_cluster = None
        single_element_cluster = 0

        for line in clusters_content:
            line = line.strip()
            if line.startswith("Cluster"):
                parts = line.split(",")
                num_of_seqs_part = parts[1].strip()
                num_of_seqs = int(num_of_seqs_part.split(":")[1].strip())
                if num_of_seqs > 1:
                    current_cluster = int(line.split()[1].replace(",", ""))
                else:
                    current_cluster = "no"
            elif line.startswith(">"):
                sequence_id = line[1:]
                cluster_mapping[sequence_id] = current_cluster

        column_name = "Predicted Cluster for h = {}".format(heights[height])
        excel_data[column_name] = excel_data["Sequence ID"].apply(lambda x: cluster_mapping.get("_" + x, "Not Found"))
        excel_data.to_excel(excel_file, index=False)
        print("DONE")