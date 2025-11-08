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

    def display_details(self):
        print(f"Sequence: {self.sequence_string} | Type: {self.sequence_type}")

    def get_alphabet(self):
        type_upper = self.sequence_type.upper()
        if type_upper == "DNA":
            return "ACGT"
        elif type_upper == "RNA":
            return "ACGU"
        elif type_upper == "PROTEIN":
            return "ACDEFGHIKLMNPQRSTVWY"
        return None

    def is_valid(self):
        alphabet_set = set(self.get_alphabet()) if self.get_alphabet() else set()
        return all(char in alphabet_set for char in self.sequence_string)

    def to_rna(self):
        if self.sequence_type.upper() != "DNA":
            return None
        rna_string = self.sequence_string.replace("T", "U")
        return BioSequence(rna_string, "RNA")

    def get_reverse_complement(self):
        if self.sequence_type.upper() != "DNA":
            return None
        base_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        reversed_complement = ''.join(base_map.get(base, base) for base in reversed(self.sequence_string))
        return BioSequence(reversed_complement, "DNA")

    def to_protein(self, offset=0):
        if self.sequence_type.upper() != "DNA":
            return None
        amino_acid_chain = ""
        sequence_length = len(self.sequence_string)
        for position in range(offset, sequence_length - 2, 3):
            triplet = self.sequence_string[position:position + 3]
            amino_acid = codon_to_amino_acid(triplet)
            amino_acid_chain += str(amino_acid) if amino_acid is not None else ""
        return BioSequence(amino_acid_chain, "PROTEIN")

class SequenceMotif:
    def __init__(self, aligned_sequences=None, weight_matrix=None, char_set=None):
        self.count_matrix = None
        if aligned_sequences:
            self.motif_length = len(aligned_sequences[0])
            self.aligned_sequences = aligned_sequences  # List of BioSequence objects
            self.char_set = aligned_sequences[0].get_alphabet()
            self._compute_counts()
            self._compute_weight_matrix()
        else:
            self.weight_matrix = weight_matrix
            self.motif_length = len(weight_matrix[0]) if weight_matrix else 0
            self.char_set = char_set

    def __len__(self):
        return self.motif_length

    def _compute_counts(self):
        num_chars = len(self.char_set)
        counts = initialize_zero_matrix(num_chars, self.motif_length)
        for sequence in self.aligned_sequences:
            for position in range(self.motif_length):
                char_index = self.char_set.index(sequence[position])
                counts[char_index][position] += 1
        self.count_matrix = counts

    def _compute_weight_matrix(self):
        if self.count_matrix is None:
            self._compute_counts()
        num_sequences = len(self.aligned_sequences)
        num_chars = len(self.char_set)
        pwm = initialize_zero_matrix(num_chars, self.motif_length)
        for i in range(num_chars):
            for j in range(self.motif_length):
                pwm[i][j] = self.count_matrix[i][j] / num_sequences
        self.weight_matrix = pwm

    def get_consensus(self):
        consensus_string = ""
        num_chars = len(self.char_set)
        for col in range(self.motif_length):
            highest_freq = -1
            best_index = -1
            for row in range(num_chars):
                if self.count_matrix[row][col] > highest_freq:
                    highest_freq = self.count_matrix[row][col]
                    best_index = row
            consensus_string += self.char_set[best_index]
        return consensus_string

    def get_threshold_consensus(self):
        consensus_string = ""
        num_sequences = len(self.aligned_sequences)
        threshold = num_sequences / 2
        num_chars = len(self.char_set)
        for col in range(self.motif_length):
            highest_freq = -1
            best_index = -1
            for row in range(num_chars):
                if self.count_matrix[row][col] > highest_freq:
                    highest_freq = self.count_matrix[row][col]
                    best_index = row
            if highest_freq > threshold:
                consensus_string += self.char_set[best_index]
            else:
                consensus_string += "-"
        return consensus_string

    def calculate_sequence_likelihood(self, target_sequence):
        likelihood = 1.0
        for pos in range(self.motif_length):
            char_index = self.char_set.index(target_sequence[pos])
            likelihood *= self.weight_matrix[char_index][pos]
        return likelihood

    def calculate_site_probabilities(self, long_sequence):
        probabilities = []
        seq_length = len(long_sequence)
        for start_pos in range(seq_length - self.motif_length + 1):
            sub_sequence = long_sequence[start_pos:start_pos + self.motif_length]
            probabilities.append(self.calculate_sequence_likelihood(sub_sequence))
        return probabilities

    def find_best_match_position(self, long_sequence):
        best_likelihood = -1.0
        best_start = -1
        seq_length = len(long_sequence)
        for start_pos in range(seq_length - self.motif_length + 1):
            sub_sequence = long_sequence[start_pos:start_pos + self.motif_length]
            current_likelihood = self.calculate_sequence_likelihood(sub_sequence)
            if current_likelihood > best_likelihood:
                best_likelihood = current_likelihood
                best_start = start_pos
        return best_start

    def align_new_motif(self, input_sequences):
        extracted_alignments = []
        for sequence in input_sequences:
            match_start = self.find_best_match_position(sequence)
            extracted_sub = sequence[match_start:match_start + self.motif_length]
            extracted_alignments.append(extracted_sub)
        return SequenceMotif(extracted_alignments)

def run_demo():
    sample1 = BioSequence("AAAGTT")
    sample2 = BioSequence("CACGTG")
    sample3 = BioSequence("TTGGGT")
    sample4 = BioSequence("GACCGT")
    sample5 = BioSequence("AACCAT")
    sample6 = BioSequence("AACCCT")
    sample7 = BioSequence("AAACCT")
    sample8 = BioSequence("GAACCT")
    all_samples = [sample1, sample2, sample3, sample4, sample5, sample6, sample7, sample8]
    motif_model = SequenceMotif(all_samples)

    print("Count Matrix:")
    display_matrix(motif_model.count_matrix)
    print("Position Weight Matrix:")
    display_matrix(motif_model.weight_matrix)
    print("Alphabet:", motif_model.char_set)
    print(motif_model.calculate_sequence_likelihood("AAACCT"))
    print(motif_model.calculate_sequence_likelihood("ATACAG"))
    print(motif_model.find_best_match_position("CTATAAACCTTACATC"))
    for seq in all_samples:
        print(seq)
    print("Consensus:", motif_model.get_consensus())
    print("Threshold Consensus:", motif_model.get_threshold_consensus())

    extended1 = BioSequence("TAAAGTTATGA")
    extended2 = BioSequence("ATGACACGTG")
    extended3 = BioSequence("TTTGGGTAT")
    updated_motif = motif_model.align_new_motif([extended1, extended2, extended3])
    display_matrix(updated_motif.count_matrix)
    print(updated_motif.find_best_match_position("AAAACT"))

if __name__ == '__main__':
    run_demo()
