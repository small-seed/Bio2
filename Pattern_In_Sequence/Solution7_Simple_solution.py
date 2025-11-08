class NaiveStringMatcher:
    @staticmethod
    def find_first_occurrence(sequence, target_pattern):
        seq_len = len(sequence)
        pat_len = len(target_pattern)
        if pat_len == 0:
            return 0
        if seq_len < pat_len:
            return -1
        
        i = 0
        found = False
        while i <= seq_len - pat_len and not found:
            j = 0
            while j < pat_len and target_pattern[j] == sequence[i + j]:
                j += 1
            if j == pat_len:
                found = True
            else:
                i += 1
        return i if found else -1
    
    @staticmethod
    def find_all_occurrences(sequence, target_pattern):
        seq_len = len(sequence)
        pat_len = len(target_pattern)
        matches = []
        if pat_len == 0:
            return list(range(seq_len + 1))
        if seq_len < pat_len:
            return []
        
        for i in range(seq_len - pat_len + 1):
            j = 0
            while j < pat_len and target_pattern[j] == sequence[i + j]:
                j += 1
            if j == pat_len:
                matches.append(i)
                i += pat_len - 1  # Skip to after the match for non-overlapping
        return matches

def interactive_test():
    input_sequence = input("Enter sequence: ").strip()
    search_pattern = input("Enter pattern: ").strip()
    
    if not search_pattern:
        print("No pattern provided.")
        return
    
    all_positions = NaiveStringMatcher.find_all_occurrences(input_sequence, search_pattern)
    
    if all_positions:
        print(f"'{search_pattern}' found in the sequence at positions: {all_positions}")
    else:
        print(f"'{search_pattern}' not found in the sequence.")

if __name__ == "__main__":
    interactive_test()
