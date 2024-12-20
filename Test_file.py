#Test
#Definieren der Nukleotide
#snakify
#from colorama import Fore, Back, Style, init
from msilib import sequence
from Bio import SeqIO
import sys


def count_nucleotides(dna_sequence_1):
    g_seq = 0
    t_seq = 0
    a_seq = 0
    c_seq = 0

    # Durchlaufen der DNA-Sequenz und Zählen der Nukleotide (addieren)
    for nucleotide in dna_sequence_1:
        if nucleotide == 'G':
            g_seq += 1
        elif nucleotide == 'T':
            t_seq += 1
        elif nucleotide == 'A':
            a_seq += 1
        elif nucleotide == 'C':
            c_seq += 1

    return g_seq, t_seq, a_seq, c_seq

def count_nucleotides(dna_sequence_2):
    g_seq_2 = 0
    t_seq_2 = 0
    a_seq_2 = 0
    c_seq_2 = 0

    for nucleotide in dna_sequence_2:
        if nucleotide == 'G':
            g_seq_2 += 1
        elif nucleotide == 'T':
            t_seq_2 += 1
        elif nucleotide == 'A':
            a_seq_2 += 1
        elif nucleotide == 'C':
            c_seq_2 += 1

    return g_seq_2, t_seq_2, a_seq_2, c_seq_2



def read_dna_file1(file_path):
    with open ("sequence1.txt", "r") as f:
        sequence = f.read().strip()
    return sequence
#öffnet die erste Datei und liest sie, lässt dabei Leerzeichen weg.
def read_dna_file2(file_path):    
    with open ("sequence2.txt", "r") as f:
        sequence = f.read().strip()
    return sequence
#öffnet die zweite Datei und liest sie, lässt dabei Leerzeichen weg.

file1 = "sequence1.txt"
file2 = "sequence2.txt"

# Beispiel DNA-Sequenz
dna_sequence_1 = read_dna_file1(file1)
dna_sequence_2 = read_dna_file2(file2)
g_seq, t_seq, a_seq, c_seq = count_nucleotides(dna_sequence_1)
g_seq_2, t_seq_2, a_seq_2, c_seq_2 = count_nucleotides(dna_sequence_2)

#vergleichen der Länge beider Sequenzen
#if len(dna_sequence_1) != len(dna_sequence_2):
    #raise ValueError("Sequenzen muessen gleich lang sein.")
#falls unterschied in Sequenz -> Mutation

if dna_sequence_1 == dna_sequence_2:
    print("Keine Mutation")
else:
    print("Mutation vorhanden") 
    #print(Fore.RED + 'Mutation vorhanden'  + Style.RESET_ALL)

def compare_dna_sequences(dna_sequence_1, dna_sequence_2):
     
    #start einer liste
    differences = []

    #Vergleich Nukleotide
    for nucleotide in range(len(dna_sequence_1)):
        if dna_sequence_1[nucleotide] != dna_sequence_2[nucleotide]:
            differences.append((nucleotide,dna_sequence_1[nucleotide], dna_sequence_2[nucleotide]))

    if differences:
        return f"Es gibt {len(differences)} Unterschiede:\n" + "\n".join([f"Position {pos}: {nuc1} -> {nuc2}" for pos, nuc1, nuc2 in differences])
    else:
        return "DNA Sequenzen sind identisch."
        
result = compare_dna_sequences(dna_sequence_1, dna_sequence_2)



print(f"{dna_sequence_1}")

print(f"G: {g_seq}")
print(f"T: {t_seq}")
print(f"A: {a_seq}")
print(f"C: {c_seq}")

print(f"{dna_sequence_2}")

print(f"G: {g_seq_2}")
print(f"T: {t_seq_2}")
print(f"A: {a_seq_2}")
print(f"C: {c_seq_2}")

print(result)