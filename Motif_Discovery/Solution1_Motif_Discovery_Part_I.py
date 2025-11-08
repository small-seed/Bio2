class Sequence:
    """
    A class to represent biological sequences (DNA, RNA, or Protein).
    Supports validation, transcription, reverse complement, and translation.
    """
    
    def __init__(self, sequence, sequence_type="DNA"):
        """
        Initialize the sequence object.
        
        :param sequence: The biological sequence as a string.
        :param sequence_type: Type of sequence ("DNA", "RNA", or "PROTEIN").
        """
        self.sequence = sequence.upper()
        self.sequence_type = sequence_type
    
    def __len__(self):
        """Return the length of the sequence."""
        return len(self.sequence)
    
    def __getitem__(self, index):
        """Support indexing to access individual characters."""
        return self.sequence[index]
    
    def __str__(self):
        """Return the string representation of the sequence."""
        return self.sequence
    
    def get_type(self):
        """Return the type of the sequence."""
        return self.sequence_type
    
    def display_info(self):
        """Print the sequence and its type."""
        print(f"Sequence: {self.sequence} | Type: {self.sequence_type}")
    
    def get_alphabet(self):
        """Return the expected alphabet for the sequence type."""
        alphabets = {
            "DNA": "ACGT",
            "RNA": "ACGU",
            "PROTEIN": "ACDEFGHIKLMNPQRSTVWY"
        }
        return alphabets.get(self.sequence_type)
    
    def is_valid(self):
        """Validate if the sequence contains only valid characters for its type."""
        alphabet = self.get_alphabet()
        if not alphabet:
            return False
        return all(char in alphabet for char in self.sequence)
    
    def transcribe(self):
        """Transcribe DNA to RNA by replacing T with U."""
        if self.sequence_type != "DNA":
            return None
        rna_seq = self.sequence.replace("T", "U")
        return Sequence(rna_seq, "RNA")
    
    def reverse_complement(self):
        """Compute the reverse complement for DNA sequences."""
        if self.sequence_type != "DNA":
            return None
        
        complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
        rev_comp = "".join(complement_map.get(base, base) for base in reversed(self.sequence))
        return Sequence(rev_comp, "DNA")
    
    def translate_to_protein(self, start_position=0):
        """
        Translate DNA sequence to protein starting from a given position.
        
        :param start_position: Starting index for translation (multiple of 3 recommended).
        :return: Sequence object of type PROTEIN.
        """
        if self.sequence_type != "DNA":
            return None
        
        protein_seq = ""
        codon_table = {
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
            "ATG": "M",
            "AAT": "N", "AAC": "N",
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
        
        for i in range(start_position, len(self.sequence) - 2, 3):
            codon = self.sequence[i:i+3]
            amino_acid = codon_table.get(codon)
            if amino_acid is None:
                protein_seq += "X"  # Use X for unknown codons
            else:
                protein_seq += amino_acid
        
        return Sequence(protein_seq, "PROTEIN")
    
    def translate_codon(self, codon):
        """
        Translate a single codon to an amino acid using the standard genetic code.
        
        :param codon: Three-letter codon string.
        :return: Amino acid letter or None if invalid.
        """
        codon_table = {
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
            "ATG": "M",
            "AAT": "N", "AAC": "N",
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
        return codon_table.get(codon.upper())
