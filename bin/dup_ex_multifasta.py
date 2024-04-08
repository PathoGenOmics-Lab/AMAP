import pandas as pd
from Bio import SeqIO
import argparse

def extraer_secuencias(fichero_tsv, fichero_fasta, fichero_salida):
    # Leer el archivo TSV
    df = pd.read_csv(fichero_tsv, sep='\t')

    # Extraer los IDs de la segunda columna y limpiarlos
    ids = []
    for row in df['Members']:
        ids.extend(row.strip("[]'").split("', '"))

    # Limpiar los IDs para extraer solo la parte relevante antes del punto y coma

    # Leer el archivo FASTA y guardar las secuencias que coinciden con los IDs
    secuencias = {record.id: record for record in SeqIO.parse(fichero_fasta, 'fasta')}
    secuencias_encontradas = [secuencias[id] for id in ids if id in secuencias]

    # Guardar las secuencias encontradas en un nuevo archivo FASTA
    SeqIO.write(secuencias_encontradas, fichero_salida, 'fasta')


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("tsv", help="tsv file")
    parser.add_argument("cds", help="cds multifasta")
    args = parser.parse_args()
    extraer_secuencias(args.tsv, args.cds, args.cds+".dup")

if __name__ == '__main__':
    main()