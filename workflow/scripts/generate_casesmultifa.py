import argparse
from Bio import SeqIO
import os
import os
import argparse
from Bio import SeqIO

class RecordProcessor:
    def __init__(self, lista_tsv_original, stay_cse1, o_file):
        self.lista_tsv_original = lista_tsv_original
        self.stay_cse1 = stay_cse1
        self.o_file = o_file

    def _extract_sequence_id(self, record):
        return record.id.split("-")[0]

    def write_to_fasta(self, record, descriptor):
        seq_name = f">{descriptor}"
        self.o_file.write(seq_name + "\n" + str(record.seq) + "\n")

    def process(self, record, case_number):
        sequence_name = self._extract_sequence_id(record)
        for entry in self.lista_tsv_original:
            if sequence_name == entry[5]:
                descriptor = self._generate_descriptor(entry, case_number)
                if self._should_write_sequence(sequence_name, case_number):
                    self.write_to_fasta(record, descriptor)

    def _generate_descriptor(self, entry, case_number):
        return f"{entry[5]};{entry[2]},{entry[3]};{entry[4]};case{case_number}"

    def _should_write_sequence(self, sequence_name, case_number):
        is_case1_and_valid = case_number == 1 and sequence_name in self.stay_cse1
        return is_case1_and_valid or case_number > 1

def read_from_file(path, condition, parsing_function):
    with open(path, "r") as infile:
        return [parsing_function(line) for line in infile if condition(line)]

def read_tsv_without_hash(path):
    return read_from_file(path, lambda line: "#" not in line, lambda line: line.strip("\n").split("\t"))

def read_first_column_without_remove(path):
    return read_from_file(path, lambda line: "remove" not in line, lambda line: line.strip("\n").split("\t")[0].split("-")[0])

def main(tsv, path2):
    input_prefix = os.path.basename(tsv).split(".")[0]
    lista_tsv_original = read_tsv_without_hash(path2 + f"{input_prefix}.tsv")
    input_prefix2 = input_prefix.split("_")[1]
    stay_cse1 = read_first_column_without_remove(f"{input_prefix2}_case1.final.tsv")

    with open(f"{input_prefix2}-complete.fa", "w") as o_file:
        processor = RecordProcessor(lista_tsv_original, stay_cse1, o_file)
        for case_number in range(1, 5):
            for record in SeqIO.parse(f"{input_prefix2}_case{case_number}.bakta.fa", "fasta"):
                processor.process(record, case_number)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Procesar archivos con prefijo variable.")
    parser.add_argument('-n', dest="name", help="Fasta name")
    parser.add_argument('-p', dest="path", help="Ruta al archivo Bakta multi-fasta.")
    parser.add_argument('-p2', dest="path2", help="Ruta al archivo Bakta multi-fasta.")
    args = parser.parse_args()
    current = os.getcwd()
    os.chdir(args.path)
    main(args.name, current + "/" + args.path2)
