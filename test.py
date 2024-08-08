#Test
#Definieren der Nukleotide
from colorama import Fore, Back, Style, init

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

# Beispiel DNA-Sequenz
dna_sequence_1 = "CGTGCATGCGTGATGATATTGCGTCGATGAC"
dna_sequence_2 = "CGTGCATGCGTGATGATATTGCGTCGATGAA"
g_seq, t_seq, a_seq, c_seq = count_nucleotides(dna_sequence_1)
g_seq_2, t_seq_2, a_seq_2, c_seq_2 = count_nucleotides(dna_sequence_2)

#vergleichen der Länge beider Sequenzen
if len(dna_sequence_1) != len(dna_sequence_2):
    raise ValueError("Sequenzen muessen gleich lang sein.")
#falls unterschied in Sequenz -> Mutation
if dna_sequence_1 == dna_sequence_2:
    print("Keine Mutation")
else:
    print("Mutation vorhanden") 
     #print(Fore.RED + "Mutation vorhanden"  + Style.RESET_ALL) text sollte rot sein, funktioniert nicht??

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