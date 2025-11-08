import re

class RestrictionEnzymeAnalyzer:
    def __init__(self, dna_sequence, enzyme_pattern):
        self.dna_sequence = dna_sequence.upper()
        self.enzyme_pattern = enzyme_pattern.upper()
    
    def _iub_to_regex(self, iub_pattern):
        iub_mapping = {
            'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
            'R': '[GA]', 'Y': '[CT]', 'M': '[AC]', 'K': '[GT]',
            'S': '[GC]', 'W': '[AT]', 'B': '[CGT]', 'D': '[AGT]',
            'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
        }
        clean_pattern = iub_pattern.replace('^', '')
        regex_pattern = ''.join(iub_mapping.get(char, char) for char in clean_pattern)
        return regex_pattern
    
    def identify_cut_sites(self):
        cut_offset = self.enzyme_pattern.find('^')
        if cut_offset == -1:
            raise ValueError("Enzyme pattern must contain '^' to indicate cut position.")
        
        regex_pat = self._iub_to_regex(self.enzyme_pattern)
        matches = re.finditer(regex_pat, self.dna_sequence)
        cut_sites = [match.start() + cut_offset for match in matches]
        return sorted(cut_sites)
    
    def generate_fragments(self, cut_sites):
        if not cut_sites:
            return [self.dna_sequence]
        
        sorted_cuts = sorted(cut_sites)
        boundaries = [0] + sorted_cuts + [len(self.dna_sequence)]
        fragments = []
        for i in range(len(boundaries) - 1):
            fragment = self.dna_sequence[boundaries[i]:boundaries[i + 1]]
            fragments.append(fragment)
        return fragments

def interactive_analysis():
    dna_input = input("Enter DNA sequence: ").strip()
    enzyme_pat = input("Enter enzyme pattern (IUB with ^ for cut): ").strip()
    
    if not dna_input or not enzyme_pat:
        print("Invalid input: Both sequence and pattern are required.")
        return
    
    analyzer = RestrictionEnzymeAnalyzer(dna_input, enzyme_pat)
    
    try:
        cut_positions = analyzer.identify_cut_sites()
        fragments = analyzer.generate_fragments(cut_positions)
        
        print(f"Cut positions: {cut_positions}")
        print(f"Resulting fragments: {fragments}")
    except ValueError as e:
        print(f"Analysis error: {e}")

if __name__ == "__main__":
    interactive_analysis()
