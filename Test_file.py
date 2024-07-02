#Hello world!
#Test
#Definieren der Nukleotide
def count_nucleotides(dna_sequence):
    g_seq = 0
    t_seq = 0
    a_seq = 0
    c_seq = 0

    # Durchlaufen der DNA-Sequenz und Zählen der Nukleotide (addieren)
    for nucleotide in dna_sequence:
        if nucleotide == 'G':
            g_seq += 1
        elif nucleotide == 'T':
            t_seq += 1
        elif nucleotide == 'A':
            a_seq += 1
        elif nucleotide == 'C':
            c_seq += 1

    return g_seq, t_seq, a_seq, c_seq

# Beispiel DNA-Sequenz
dna_sequence = "AGCTTAGCCTAGGCTAAGATAGTGCTGCTCGATGTAGTCGGCGCGAGAGAGCTTCTCATCGTGCATGCGTGATGATATTGCGTCGATGAC"
g_seq, t_seq, a_seq, c_seq = count_nucleotides(dna_sequence)

print(f"{dna_sequence}")
print(f"G: {g_seq}")
print(f"T: {t_seq}")
print(f"A: {a_seq}")
print(f"C: {c_seq}")