import sys
from Bio import SeqIO

fasta_file = "gene.fna" #Pfad zur FASTA Datei


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

    total_nucleotides = g_seq + t_seq + a_seq + c_seq
    
    return g_seq, t_seq, a_seq, c_seq, total_nucleotides




# Die FASTA-Datei öffnen und die Sequenzen lesen
with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        # Zähle die Nukleotide in der aktuellen Sequenz
        g_seq, t_seq, a_seq, c_seq, total_nucleotides = count_nucleotides(str(record.seq))

        # Ergebnisse für die aktuelle Sequenz ausgeben
        print(f"Sequenz-ID: {record.id}")
        print(f"Anzahl von A: {a_seq}")
        print(f"Anzahl von T: {t_seq}")
        print(f"Anzahl von C: {c_seq}")
        print(f"Anzahl von G: {g_seq}")
        print(f"Total: {total_nucleotides} ")
        print("----------")


#Damit die Sequenzen korrekt miteinander verglichen werden muss in Abschnitten gearbeitet werden. Eine Insertion oder Delition kann zu einer kompletten Verschiebung des Codes führen. 
#Die Frage ist jetzt wie soll die Ausgabe sein, dass sie 1. nicht zu lang ist und 2. dass man versteht wo die Unterschiede liege und was diese aussagen

