import sys
from Bio import SeqIO
#test 1
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

def compare_dna_sequences(dna_sequence_1, dna_sequence_2):
    min_length = min(len(dna_sequence_1), len (dna_sequence_2))
    #sucht sich die kürzere Sequenz

    #vergleicht Sequenzen bis zu Länge der kürzeren Sequenz
    #i erhöht sich bei jedem Durchgang um 1
    for i in range(min_length):
        if dna_sequence_1[i] != dna_sequence_2[i]:
            return f"Erster Unterschied ist an Position {i + 1}: {dna_sequence_1[i]} -> {dna_sequence_2[i]}"
    
    return "Die beiden Sequenzen stimmen bis zur Länge der kürzeren Sequenz überein."

sequences = []

# Die FASTA-Datei öffnen und die Sequenzen lesen
with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
# Sequenzen sammeln
        #str(record.seq) wandelt Seq-Objekt in normalen String
        #append fügt das Ergebnis von str(record.seq) an Ende der Liste sequences hinzu
        sequences.append(str(record.seq))

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

#überprüfen ob genau 2 sequenzen in Datei
if len(sequences) == 2:
    #liste beginnt bei 0
    dna_sequence_1 = sequences[0]
    dna_sequence_2 = sequences[1]
    #Zeige den ersten Unterschied an
    result = compare_dna_sequences(dna_sequence_1, dna_sequence_2)
    print(result)
else:
    print("Die Datei enthaelt nicht genau zwei Sequenzen.")

#Notizen:
#Damit die Sequenzen korrekt miteinander verglichen werden muss in Abschnitten gearbeitet werden. Eine Insertion oder Delition kann zu einer kompletten Verschiebung des Codes führen. 
#Die Frage ist jetzt wie soll die Ausgabe sein, dass sie 1. nicht zu lang ist und 2. dass man versteht wo die Unterschiede liege und was diese aussagen
#Eine Idee wäre es nicht die Unterscheide, sondern die Gemeinsamkeiten zu finden, um zu wissen wo die Mutation ähnlich sind und was sie ausmachen



