import sys
from Bio import SeqIO
import pandas as pd
#Pandas is a Python library. Pandas is used to analyze data. (W3Schools)
#Pandas wird verwendet um nicht mit Listen zu arbeiten sondern mit Tabellen.

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
   
    """
    diese Funktion vergleicht gleich (part_length) lange Teile von 2 Sequenzen, 
    angefagen von den jeweils angegeben start_index_1 bzw. start_index_2. Zurückgegeben wird die Anzahl Unterschiede.
    """

    mismatch_count = 0
    for k in range(part_length):
        if start_index_1 + k < len(dna_sequence_1) and start_index_2 + k < len(dna_sequence_2):
            if dna_sequence_1[start_index_1 + k] != dna_sequence_2[start_index_2 + k]:
                mismatch_count += 1
    return mismatch_count



def compare_dna_sequences(dna_sequence_1, dna_sequence_2):
    i, j = 0, 0 #Braucht zwei Variablen wegen Deletion
    differences = []
    #dictionary zur Speicherung pathogener Mutationen
    pathogen_mutations = {
        76190: ("Substitution", "G", "A"),
        121388: ("Deletion", "C", "_"),
    }

    #i und j erhöht sich bei jedem Durchgang um 1
    while i < len(dna_sequence_1) and j < len(dna_sequence_2): 
    #while schleife, solange es in der Grenze beider Sequenzen liegt, wenn Länge überschritten Schliefe hört auf
        if dna_sequence_1[i] != dna_sequence_2[j]:
            #überprüfen auf deletion

            #listen für Mismatch-Zählungen
            mismatch_count_i = [compare_dna_sequence_parts(dna_sequence_1, dna_sequence_2, i + k, j, 10) for k in range(6)]
            mismatch_count_j = [compare_dna_sequence_parts(dna_sequence_1, dna_sequence_2, i, j + k, 10) for k in range(6)]

            #suche nach dem ersten k, bei dem mismatch_count_i oder mismatch_count_j 0 ist
            for k in range(6):
                if mismatch_count_i[k] == 0 or mismatch_count_j[k] == 0:
                    if mismatch_count_i[k] == 0:
                        #füge einen Eintrag in differnces Liste hinzu, der Deletion beschreibt
                        differences.append({
                            "Art": "Deletion", #Mutationstyp
                            "Position": i + 1, #Position der Mutation
                            "Original": dna_sequence_1[i], #Position ohne Mutation
                            "Mutation": "_", #Deletion ist _
                            "Betroffene Nukleotide": k, #wie viele Nukleotide betroffen sind
                            "Auswirkung": "pathogen" if pathogen_mutations.get(i + 1) == ("Deletion", dna_sequence_1[i], "_") else "" #schaut ob diese Mutation als pathgen definiert ist in pathogen Dictionary
                        })
                        i += k
                    if mismatch_count_j[k] == 0:
                        j += k
                    break
            else:
                differences.append({
                    #füge Eintrag in differences Liste hinzu, der Substitution beschreibt
                    #gleiches vorgehen wie bei Deletion
                    "Art": "Substitution",
                    "Position": i + 1,
                    "Original": dna_sequence_1[i],
                    "Mutation": dna_sequence_2[j],
                    "Betroffene Nukleotide": 1,
                    "Auswirkung": "pathogen" if pathogen_mutations.get(i + 1) == ("Substitution", dna_sequence_1[i], dna_sequence_2[j]) else "" #schaut ob diese Mutation als pathgen definiert ist in pathogen Dictionary
                })
                i += 1
                j += 1
        else:
            #wenn Nukleotide übereinstimmen, einfach weiter
            i += 1
            j += 1

    while i < len(dna_sequence_1):
        #wenn Wildtyp länger ist dann sind am ende noch Mutationen die hier berücksichtigt werden
        #nur notwendig wenn letztes Nukleotid gelöscht wurde
        differences.append({
            "Art": "Deletion",
            "Position": i + 1,
            "Original": dna_sequence_1[i],
            "Mutation": "_",
            "Betroffene Nukleotide": 1, 
            "Auswirkung": "" #keine spezifischen pathogenen Mutationen definiert
        })
        i += 1

    return differences
    

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

differences = compare_dna_sequences(sequences_1[0], sequences_2[0])
df = pd.DataFrame(differences) #erstellt eine Tabelle aus Liste von Dictionaries differnces

# gibt Tabelle formartiert aus, Spalten sind linksbündig, Index wird weggelassen
print(df.to_string(index=False, justify='left'))