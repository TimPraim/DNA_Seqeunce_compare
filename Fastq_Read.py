from Bio import SeqIO

file_path = "C:/Users/clarb/OneDrive - Departement Bildung, Kultur und Sport Kanton Aargau-7988362-Kantonsschule Zofingen/Matura/DNA Sequenzen/SRR390728.fastq/SRR390728.fastq"


for record in  SeqIO.parse(file_path,"fastq"):
    print(f"ID: {record.id}")


