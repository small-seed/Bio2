import re

class ProteinMotifRegexDetector:    
    def __init__(self, amino_sequence, prosite_pattern):
        self.amino_sequence = amino_sequence.upper()
        self.prosite_pattern = prosite_pattern
    
    def detect_zinc_finger_ring(self):
        zinc_regex = r'C\.[H][LIVMFY]C\.{2}C[LIVMYA]'
        match = re.search(zinc_regex, self.amino_sequence)
        return match.start() if match else -1
    
    def _convert_prosite_to_regex(self, pattern):
        regex_pat = re.sub(r'−|-', '', pattern)
        # Replace x (any single) with .
        regex_pat = regex_pat.replace('x', '.')
        # Handle repeats: x(n) -> .{n}, but since x replaced, (n) -> {n}
        regex_pat = re.sub(r'\((\d+)\)', r'{\1}', regex_pat)
        
        return regex_pat
    
    def detect_prosite_motif(self):
        regex_pattern = self._convert_prosite_to_regex(self.prosite_pattern)
        match = re.search(regex_pattern, self.amino_sequence)
        return match.start() if match else -1

def interactive_test():
    protein_seq = input("Enter protein sequence: ").strip()
    prosite_pat = input("Enter ProSite pattern: ").strip()
    
    detector = ProteinMotifRegexDetector(protein_seq, prosite_pat)
    
    # Detect Zinc finger
    zinc_pos = detector.detect_zinc_finger_ring()
    print(f"Zinc finger RING motif found at position: {zinc_pos}")
    
    # Detect general ProSite
    prosite_pos = detector.detect_prosite_motif()
    print(f"ProSite motif found at position: {prosite_pos}")

if __name__ == "__main__":
    interactive_test()
