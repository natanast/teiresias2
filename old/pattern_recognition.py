from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

# Read the DNA sequences from a FASTA file
sequences = list(SeqIO.parse("sequences.fasta", "fasta"))

# Create a MultipleSeqAlignment object
alignment = MultipleSeqAlignment(sequences)

# Use Clustalw to perform the alignment
clustalw_cline = ClustalwCommandline("clustalw", infile="sequences.fasta")
stdout, stderr = clustalw_cline()

# Save the alignment to a FASTA file
SeqIO.write(alignment, "aligned_sequences.fasta", "fasta")

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# Read the DNA sequences from the text file
with open("sequences.txt", "r") as handle:
    sequences = handle.read().splitlines()

# Create a MultipleSeqAlignment object with the DNA sequences
alignment = MultipleSeqAlignment(sequences)

# Use the align function from the MultipleSeqAlignment object to perform the alignment
aligned_alignment = alignment.align()

# Print the aligned alignment
print(aligned_alignment)
