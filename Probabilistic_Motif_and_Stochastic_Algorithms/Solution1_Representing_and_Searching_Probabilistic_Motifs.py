def initialize_zero_matrix(row_count, col_count):
    return [[0 for _ in range(col_count)] for _ in range(row_count)]

def display_matrix(matrix):
    for row in matrix:
        print(row)

def codon_to_amino_acid(codon):
    code_table = {
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TGT": "C", "TGC": "C",
        "GAT": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "TTT": "F", "TTC": "F",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAT": "H", "CAC": "H",
        "ATA": "I", "ATT": "I", "ATC": "I",
        "AAA": "K", "AAG": "K",
        "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATG": "M", "AAT": "N", "AAC": "N",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TGG": "W",
        "TAT": "Y", "TAC": "Y",
        "TAA": "_", "TAG": "_", "TGA": "_"
    }
    return code_table.get(codon, None)

class BioSequence:
    def __init__(self, sequence_string, sequence_type="DNA"):
        self.sequence_string = sequence_string.upper()
        self.sequence_type = sequence_type

    def __len__(self):
        return len(self.sequence_string)

    def __getitem__(self, index_or_slice):
        return self.sequence_string[index_or_slice]

    def __str__(self):
        return self.sequence_string

    def sequence_type(self):
        return self.sequence_type

    def display
