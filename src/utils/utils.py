from functools import wraps
from time import time


def read_file(file_path: str):
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


def timeit(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time()
        result = func(*args, **kwargs)
        end = time()
        print(f"{func.__name__} took {end - start} seconds")
        return result

    return wrapper
