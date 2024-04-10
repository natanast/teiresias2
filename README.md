# Hierarchical Clustering for Protein Sequences - MSc Thesis

This Python code performs hierarchical clustering on input protein sequences, allowing users to analyze the similarities between sequences and group them into clusters based on those similarities. The resulting clusters are then written to an Excel file for further analysis.

## Installation  

1. First create a virtual enviroment 
    ```bash
    python -m venv venv 
    venv\\Scripts\\activate (for Windows) 
    source venv/bin/activate (for Linux/Mac) 
    ```

2. Install the dependencies
    ```bash
    pip install -r requirements.txt
    ```

## Usage

```bash
python main.py -i <input> [-m <multiplier>] [-w <weight>] [-p <penalty>] [-t <threshold>] [-c <cuts>] -o <output>
```
The script `main.py` accepts the following command-line arguments:

- `-i, --input`: Path to the input file containing amino acid (AA) sequences.
- `-m, --multiplier`: Multiplier for similarity calculation after a specified threshold amino acid position.
- `-w, --weight`: Added weight for concurrent matches.
- `-p, --penalty`: Penalty for mismatch.
- `-t, --threshold`: Amino acid position number after which the match will be multiplied with the multiplier argument.
- `-c, --cuts`: Number of cuts in the hierarchical tree (4, 8, or 12).
- `-o, --output`: Path to the Excel file where the calculated clusters will be written as new columns.

## Arguments Explanation

- `input`: This argument specifies the path to the input file containing AA sequences. It's mandatory to provide this argument.
- `multiplier`: Specifies the multiplier for similarity calculation after a certain threshold amino acid position. Default value is 2.0.
- `weight`: Added weight for concurrent matches. Default value is 0.1.
- `penalty`: Penalty for mismatch. Default value is 0.0.
- `threshold`: Amino acid position number after which the match will be multiplied with the multiplier argument. Default value is 105.
- `cuts`: Number of cuts in the hierarchical tree (Can take values 4, 8, or 12). Default value is 4.
- `output`: Specifies the path to the Excel file where the calculated clusters will be written as new columns.

## Input File Format

The input file is expected to have a .fa format and contains AA sequences, each sequence on a separate line.

## Output

The script generates a distance matrix based on the input sequences and saves it as `UPGMA_Input.txt`. It then creates a hierarchical tree and divides it into clusters based on the specified number of cuts. Then, it creates a file for each tree cut (on master-thesis-main/results/clusters/) with its corresponding calculated clusters and statistics like inter and outer cluster similarity for each cluster. Finally, it writes the calculated clusters to the specified Excel file for every sequence.

## Example run

```bash
python main.py -i C:/Users/30694/Desktop/Nina/master-thesis-main/data/tiny_test.fa -m 2.0 -c 4 -p 0.0 -w 0.1 -o C:/Users/30694/Desktop/Nina/master-thesis-main/data/CLL-DB-data-aligned_14_seqs.xlsx
```

The cmd output should look like this: 

```bash
Total number of sequences: 14
Calculating distance: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 14/14 [00:00<00:00, 2275.89it/s] 
distance_matrix took 0.019102096557617188 seconds
Distance matrix saved

Creating tree...
get_tree took 0.04251408576965332 seconds
Tree created


Tree length is:  0.2942935526371002
Tree will be cut on the following heights:  [0.12609303928911686, 0.20445840433239937, 0.24568521697074175, 0.2750881714746356]
```