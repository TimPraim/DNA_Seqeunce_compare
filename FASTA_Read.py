import sys
from Bio import SeqIO
fasta_file1 = "BRCA1 normal.fna" #Pfad zur FASTA Datei
fasta_file2 = "BRCA1 mutation simulation.fna"


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

def compare_dna_sequence_parts(dna_sequence_1, dna_sequence_2, start_index_1, start_index_2, part_length):
    #diese Funktion vergleicht gleich (part_length) lange Teile von 2 Sequenzen, angefagen von den jeweils angegeben start_index_1 bzw. start_index_2. Zurückgegeben wird die Anzahl Unterschiede.
    mismatch_count = 0
    for k in range(part_length):
        if start_index_1 + k < len(dna_sequence_1) and start_index_2 + k < len(dna_sequence_2):
            if dna_sequence_1[start_index_1 + k] != dna_sequence_2[start_index_2 + k]:
                mismatch_count += 1
    return mismatch_count

def compare_next_nucleotides(dna_sequence_1, dna_sequence_2, i, j):
    #vergleicht die nächsten nukleotide ab position i in der ersten sequenz mit den nächsten nukleotide von j

       #berechne die Endindizes für den Vergleich
    end_index_1 = i + 20

    #schleife, um j so lange zu erhöhen, bis die nächsten 10 Nukleotide übereinstimmen
    while j < len(dna_sequence_2):
        end_index_2 = j + 20

        #überprüfe, ob die Indizes innerhalb der Grenzen liegen
        if end_index_1 <= len(dna_sequence_1) and end_index_2 <= len(dna_sequence_2):
        #vergleiche die nächsten 10 Nukleotide
            if dna_sequence_1[i:end_index_1] == dna_sequence_2[j:end_index_2]:
                return i, j + 1  #wenn sie übereinstimmen, j um 1 erhöhen
            else:
                #wenn sie nicht übereinstimmen, j um 2 erhöhen
                j += 2
        else:
            #falls die Grenzen überschritten werden, erhöhen wir einfach j um 2, um Fehler zu vermeiden
            j += 2

    return i, j  # Rückgabe der neuen Indizes

def compare_dna_sequences(dna_sequence_1, dna_sequence_2):
    i, j = 0, 0 #Braucht zwei Variablen wegen Deletion
    differences = []
    #i und j erhöht sich bei jedem Durchgang um 1
    while i < len(dna_sequence_1) and j < len(dna_sequence_2): 
    #while schleife, solange es in der Grenze beider Sequenzen liegt, wenn Länge überschritten Schliefe hört auf
        if dna_sequence_1[i] != dna_sequence_2[j]:
            mismatch_count = compare_dna_sequence_parts(dna_sequence_1, dna_sequence_2, i, j, 6)
            #überprüfen auf deletion

            #arrays für Mismatch-Zählungen
            mismatch_count_i = [compare_dna_sequence_parts(dna_sequence_1, dna_sequence_2, i + k, j, 10) for k in range(6)]
            mismatch_count_j = [compare_dna_sequence_parts(dna_sequence_1, dna_sequence_2, i, j + k, 10) for k in range(6)]

            #suche nach dem ersten k, bei dem mismatch_count_i oder mismatch_count_j 0 ist
            for k in range(6):
                if mismatch_count_i[k] == 0 or mismatch_count_j[k] == 0:
                    if mismatch_count_i[k] == 0:
                        #verschiebe i um k, wenn i verschoben werden kann
                        differences.append(f"Deletion an Position {i + 1}: {dna_sequence_1[i]} -> _")  #melde die Deletion
                        i += k  #i wird aktualisiert
                    if mismatch_count_j[k] == 0:
                        #verschiebe j um k, wenn j verschoben werden kann
                        j += k
                    break
            else:
                #wenn keine 0 in den mismatch_counts gefunden wurden, es könnte eine Substitution sein
                differences.append(f"Substitution an Position {i + 1}: {dna_sequence_1[i]} -> {dna_sequence_2[j]}")
                i += 1
                j += 1
        else:
            #wenn Nukleotide übereinstimmen, einfach weiter
            i += 1
            j += 1

    #überprüfen, ob noch Unterschiede in der Referenzsequenz vorhanden sind
    while i < len(dna_sequence_1):
        differences.append(f"Zusätzliche Base in der Referenzsequenz an Position {i + 1}: {dna_sequence_1[i]} fehlt in der mutierten Sequenz.")
        i += 1

    if differences:
        return "\n".join(differences)
    else:
        return "Keine Unterschiede gefunden."
    

sequences_1 = []
sequences_2 = []

# Die FASTA-Datei öffnen und die Sequenzen lesen
with open(fasta_file1, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
# Sequenzen sammeln
        #str(record.seq) wandelt Seq-Objekt in normalen String
        #append fügt das Ergebnis von str(record.seq) an Ende der Liste sequences hinzu
        sequences_1.append(str(record.seq))

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

# Die FASTA-Datei öffnen und die Sequenzen lesen
with open(fasta_file2, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
# Sequenzen sammeln
        #str(record.seq) wandelt Seq-Objekt in normalen String
        #append fügt das Ergebnis von str(record.seq) an Ende der Liste sequences hinzu
        sequences_2.append(str(record.seq))

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

# Vergleiche die beiden Sequenzen
result = compare_dna_sequences(sequences_1[0], sequences_2[0])
print(result)
