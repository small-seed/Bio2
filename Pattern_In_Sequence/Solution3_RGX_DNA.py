import re

class DNAPatternMatcher:    
    def __init__(self, dna_sequence, target_pattern):
        self.dna_sequence = dna_sequence
        self.target_pattern = target_pattern
    
    def locate_first_match(self, sequence, pattern):
        match = re.search(pattern, sequence)
        return match.start() if match else -1
    
    def locate_all_matches(self, sequence, pattern):
        matches = re.finditer(pattern, sequence)
        return [match.start() for match in matches]

def demonstrate_usage():
    dna_seq = input("Enter DNA sequence: ").strip()
    regex_pat = input("Enter pattern (regex): ").strip()
    
    matcher = DNAPatternMatcher(dna_seq, regex_pat)
    
    # Find first occurrence
    first_pos = matcher.locate_first_match()
    if first_pos >= 0:
        print(f"First match found at position: {first_pos}")
    else:
        print("No match found")
    
    # Find all occurrences
    all_positions = matcher.locate_all_matches()
    if all_positions:
        print(f"All matches found at positions: {all_positions}")
    else:
        print("No matches found")

if __name__ == "__main__":
    demonstrate_usage()
