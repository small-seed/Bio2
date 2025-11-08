import re

class SequenceRegexAnalyzer:    
    def __init__(self, input_sequence, allowed_chars="ACGT"):
        self.input_sequence = input_sequence.upper()
        self.allowed_chars = allowed_chars
    
    def validate_sequence_composition(self):
        invalid_pattern = re.escape(self.allowed_chars)
        if re.search(f'[^{invalid_pattern}]', self.input_sequence):
            return False
        return True
    
    def _get_codon_table(self):
        return {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        }
    
    def convert_to_protein(self):
        codon_mapping = self._get_codon_table()
        if len(self.input_sequence) % 3 != 0:
            raise ValueError("Sequence length must be a multiple of 3 for translation.")
        
        protein_chain = ""
        for start_idx in range(0, len(self.input_sequence), 3):
            codon = self.input_sequence[start_idx:start_idx + 3]
            if codon in codon_mapping:
                protein_chain += codon_mapping[codon]
            else:
                raise ValueError(f"Invalid codon encountered: {codon}")
        return protein_chain
    
    def identify_longest_protein_segment(self):
        try:
            protein = self.convert_to_protein()
        except ValueError:
            return ""
        
        matches = re.finditer(r'M[^_]*_', protein)
        max_length = 0
        longest_segment = ""
        
        for match_obj in matches:
            segment_start, segment_end = match_obj.span()
            segment_length = segment_end - segment_start
            if segment_length > max_length:
                max_length = segment_length
                longest_segment = match_obj.group()
        
        return longest_segment

def interactive_demo():
    dna_input = input("Enter the sequence: ").strip()
    regex_input = input("Enter pattern (as a regular expression; used for validation chars): ").strip()
    
    analyzer = SequenceRegexAnalyzer(dna_input, regex_input)
    
    # Perform translation
    try:
        translated_protein = analyzer.convert_to_protein()
        print(f"Translated protein: {translated_protein}")
    except ValueError as e:
        print(f"Translation error: {e}")
    
    # Find longest protein segment
    longest = analyzer.identify_longest_protein_segment()
    print(f"Longest protein segment: {longest}")

if __name__ == "__main__":
    interactive_demo()
