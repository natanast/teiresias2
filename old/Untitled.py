from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator

# Load the sequences from a FASTA file
sequences = list(SeqIO.parse("output.fa", "fasta"))
print(sequences)

# # Align the sequences using Clustal Omega
# clustalomega_cline = ClustalOmegaCommandline(infile="output.fa", outfile="alignment.fasta", verbose=True, auto=True)
# clustalomega_cline()

# # Load the alignment from the output file
# alignment = AlignIO.read("alignment.fasta", "fasta")

# # Calculate the pairwise distances between all sequences
calculator = DistanceCalculator('identity')
distances = calculator.get_distance(sequences)

# # Print the pairwise distances
for i, name_i in enumerate(sequences):
     for j, name_j in enumerate(sequences):
         if i < j:
             print(f"{name_i.id} vs {name_j.id}: {distances[i][j]:.4f}")


