from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from collections import Counter
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW
from Bio import SearchIO


def parse_seq_fasta(file_path):
    record = SeqIO.read(file_path, "fasta")
    return record.seq


def get_GC_content(sequence):
    return gc_fraction(sequence) * 100


def get_nucleotides_count(sequence):
    return {
        "total": len(sequence),
        "A": sequence.count("A"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
        "T": sequence.count("T"),
    }


def get_frequency_count(sequence):
    c_nucl = get_nucleotides_count(sequence)
    return {
        "A": c_nucl["A"] / c_nucl["total"] * 100,
        "G": c_nucl["G"] / c_nucl["total"] * 100,
        "C": c_nucl["C"] / c_nucl["total"] * 100,
        "T": c_nucl["T"] / c_nucl["total"] * 100,
    }


def draw_nucleotides_frequency_visualisation(sequence):
    frequency = get_frequency_count(sequence)
    plt.title("SARS-Cov2 Nucleotides frequency")
    plt.pie(frequency.values(), labels=frequency.keys(), autopct='%1.1f%%')
    plt.show()


def transcribe(dna_sequence):
    return dna_sequence.transcribe()


def translate(rna_sequence):
    return pad_seq(rna_sequence).translate()


# чтобы не было ворнинга, что число нуклеотидов в последовательности
# не кратно трём
def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """
    r = len(sequence) % 3
    return sequence if r == 0 else sequence + Seq('N' * (3 - r))


def get_codons(sequence):
    codons_generator = (sequence[n:n + 3] for n in range(0, len(sequence), 3))
    return list(codons_generator)


def draw_most_common_codons_frequency(sequence):
    most_common_codons = get_most_common_entries_frequency(
        get_codons(sequence))
    plt.title('Most common codons frequency')
    plt.bar(range(len(most_common_codons)),
            most_common_codons.values(), align='center')
    plt.xticks(range(len(most_common_codons)),
               most_common_codons.keys())
    plt.ylabel("%")
    plt.xlabel("Codon")
    plt.show()


def get_amino_acids(sequence):
    amino_acids = translate(transcribe(sequence)).replace('*', '')
    return amino_acids


def get_most_common_entries_frequency(entries_list, most_common_int=10):
    most_common_entries = dict(Counter(entries_list).most_common(
        most_common_int))
    entries_frequency = {}
    for key in most_common_entries.keys():
        entries_frequency[key] = (most_common_entries[key] / len(
            entries_list)) * 100
    return entries_frequency


def draw_most_common_amino_acids_plot(sequence):
    most_common_amino_acids = get_most_common_entries_frequency(
        get_amino_acids(sequence))
    plt.title('Most common amino acids frequency')
    plt.bar(range(len(most_common_amino_acids)),
            most_common_amino_acids.values(), align='center')
    plt.xticks(range(len(most_common_amino_acids)),
               most_common_amino_acids.keys())
    plt.ylabel("%")
    plt.xlabel("Amino acid")
    plt.show()


def get_blast_info(file_name):
    seq_protein = translate(transcribe(parse_seq_fasta(file_name)))
    result_handle = NCBIWWW.qblast("blastp", "pdb", seq_protein)
    records = SearchIO.read(result_handle, 'blast-xml')
    polymerases = []
    endoribonucleases = []
    for record in records:
        desc = record.description
        if 'polymerase' in desc:
            polymerases.append(record)
        elif 'Uridylate-specific endoribonuclease' in desc:
            endoribonucleases.append(record)
        else:
            print(record.id, desc)

    print('proteins count:', len(records))
    print('polymerases count:', len(polymerases))
    print('Uridylate-specific endoribonucleases count:', len(endoribonucleases))


get_blast_info('covid_sequence.fasta')

# draw_nucleotides_frequency_visualisation(
#     parse_seq_fasta('covid_sequence.fasta'))
# draw_most_common_amino_acids_plot(parse_seq_fasta('covid_sequence.fasta'))
# draw_most_common_codons_frequency(parse_seq_fasta('covid_sequence.fasta'))
