import argparse
from Bio import SeqIO
import os

def main(tsv, path2):
    input_prefix = tsv
    
    lista_tsv_original = []
    with open(path2+f"{input_prefix}.tsv", "r") as infile:
        for line in infile:
            if "#" not in line:
                l = line.strip("\n").split("\t")
                lista_tsv_original.append(l)

    stay_cse1 = []
    with open(f"{input_prefix}_case1.final.tsv", "r") as infile:
        for line in infile:
            if "remove" not in line:
                l = line.strip("\n").split("\t")[0].split("-")[0]
                stay_cse1.append(l)

    with open(f"{input_prefix}-complete.fa", "w") as o_file:
        for case_number in range(1, 5):
            for registro in SeqIO.parse(f"{input_prefix}_case{case_number}.bakta.fa", "fasta"):
                nombre_secuencia = registro.id
                secuencia = str(registro.seq)

                for elemento in lista_tsv_original:
                    e = nombre_secuencia.split("-")[0]
                    if e == elemento[5]:
                        if case_number == 1 and e in stay_cse1:
                            seq_name = f">{elemento[5]};{elemento[2]},{elemento[3]};{elemento[4]};case{case_number}"
                            seq = seq_name + "\n" + secuencia + "\n"
                            o_file.write(seq)
                        elif case_number > 1:
                            seq_name = f">{elemento[5]};{elemento[2]},{elemento[3]};{elemento[4]};case{case_number}"
                            seq = seq_name + "\n" + secuencia + "\n"
                            o_file.write(seq)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Procesar archivos con prefijo variable.")
    parser.add_argument('-n', dest="name", help="Fasta name")
    parser.add_argument('-p', dest="path", help="Ruta al archivo Bakta multi-fasta.")
    parser.add_argument('-p2', dest="path2", help="Ruta al archivo Bakta multi-fasta.")
    args = parser.parse_args()
    current = os.getcwd()
    os.chdir(args.path)
    main(args.name, current+"/"+args.path2)
