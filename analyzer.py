from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from collections import Counter
import argparse
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW
from Bio import SearchIO


def parse_seq_fasta(file_path):
    record = SeqIO.read(file_path, "fasta")
    return record.seq


def get_GC_content(sequence):
    return gc_fraction(sequence) * 100


def count_nucleotides(sequence):
    return {
        "total": len(sequence),
        "A": sequence.count("A"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
        "T": sequence.count("T"),
    }


def frequency_nucleotides(sequence):
    c_nucl = count_nucleotides(sequence)
    return {
        "A": c_nucl["A"] / c_nucl["total"] * 100,
        "G": c_nucl["G"] / c_nucl["total"] * 100,
        "C": c_nucl["C"] / c_nucl["total"] * 100,
        "T": c_nucl["T"] / c_nucl["total"] * 100,
    }


def transcribe(dna_sequence):
    return dna_sequence.transcribe()


def translate(rna_sequence):
    return rna_sequence.translate()


# чтобы не было ворнинга, что число нуклеотидов в последовательности
# не кратно трём
# def pad_seq(sequence):
#     """ Pad sequence to multiple of 3 with N """
#     r = len(sequence) % 3
#     return sequence if r == 0 else sequence + Seq('N' * (3 - r))


def count_codons(sequence):
    codons_generator = (sequence[n:n + 3] for n in range(0, len(sequence), 3))
    codons = Counter(codons_generator).most_common()
    return codons


def most_common_amino_acids(sequence):
    amino_acids = translate(transcribe(sequence)).replace('*', '')
    return Counter(amino_acids).most_common(10)


def get_blast_info(file_name):
# pdb|8G6R|A Chain A, nsp12 [Porcine epidemic diarrhea virus]
    seq_protein = translate(transcribe(parse_seq_fasta(file_name)))
    result_handle = NCBIWWW.qblast("blastp", "pdb", seq_protein)
    records = SearchIO.read(result_handle, 'blast-xml')
    polymerases = []
    helicases = []
    for record in records:
        desc = record.description
        if 'polymerase' in desc:
            polymerases.append(record.id)
        elif 'helicase' in desc or 'Helicase' in desc:
            helicases.append(record.id)
        else:
            print(record.id, desc)

    print('polymerases count', len(polymerases))
    print('helicases count:', len(helicases))


def main():
    argparser = argparse.ArgumentParser()
    proc_arg_group = argparser.add_mutually_exclusive_group()
    proc_arg_group.add_argument('-trb', '--transcribe', nargs=1, type=str,
                                metavar='sequence_file', dest='transcribe',
                                required=False, help='Транскрипция '
                                                     'последовательности из '
                                                     'файла')
    a_group = argparser.add_mutually_exclusive_group()
    a_group.add_argument('sequence_file', type=str, action='store', nargs='?',
                         help='Путь до файла с геномными данными '
                              'в формате .fasta')
    args = argparser.parse_args()

    if args.transcribe:
        print(transcribe(parse_seq_fasta(args.transcribe[0])))
    elif args.sequence_file:
        print('Здесь будет полный анализ')
    else:
        print('Такой команды нет\n'
              'Справка -h | --help')


# if __name__ == '__main__':
#     main()

get_blast_info("covid_sequence.fasta")
