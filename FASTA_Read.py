import sys
from Bio import SeqIO

fasta_file = "gene.fna"

# Die FASTA-Datei öffnen und die Sequenz lesen
with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        # Die ersten 10 Nukleotide aus der Sequenz ausgeben
        print(record.seq[:10])
        break  # Nur den ersten Eintrag verarbeiten