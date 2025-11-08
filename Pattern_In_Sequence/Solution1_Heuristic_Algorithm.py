import os
import re

class BoyerMooreSearcher:
    def __init__(self, char_set, target_pattern):
        self.char_set = char_set
        self.target_pattern = target_pattern
        self._m = len(target_pattern)
        self._precompute_tables()
    
    def _precompute_tables(self):
        self._build_bad_char_table()
        self._build_good_suffix_table()
    
    def _build_bad_char_table(self):
        self.bad_char = {}
        for ch in self.char_set:
            self.bad_char[ch] = -1
        for pos in range(self._m):
            self.bad_char[self.target_pattern[pos]] = pos
    
    def _build_good_suffix_table(self):
        m = self._m
        self.border_pos = [0] * (m + 1)
        self.shift_val = [0] * (m + 1)
        
        # Initialize for the full pattern
        i = m
        j = m + 1
        self.border_pos[i] = j
        
        # Compute borders for strong good suffix
        while i > 0:
            while j <= m and self.target_pattern[i - 1] != self.target_pattern[j - 1]:
                if self.shift_val[j] == 0:
                    self.shift_val[j] = j - i
                j = self.border_pos[j]
            i -= 1
            j -= 1
            self.border_pos[i] = j
        
        # Handle case 2: fill remaining shifts using borders
        j = self.border_pos[0]
        for i in range(m + 1):
            if self.shift_val[i] == 0:
                self.shift_val[i] = j
            if i == j:
                j = self.border_pos[j]
    
    def _compute_bad_char_shift(self, mismatch_pos, mismatch_char):
        last_occ = self.bad_char.get(mismatch_char, -1)
        if last_occ > mismatch_pos:
            return mismatch_pos + 1
        else:
            return mismatch_pos - last_occ
    
    def find_occurrences(self, input_text):
        results = []
        n = len(input_text)
        i = 0
        while i <= n - self._m:
            j = self._m - 1
            # Compare from the end
            while j >= 0 and self.target_pattern[j] == input_text[i + j]:
                j -= 1
            if j < 0:
                # Full match found
                results.append(i)
                i += self.shift_val[0]
            else:
                # Mismatch: apply max of good suffix and bad char shifts
                mismatch_ch = input_text[i + j]
                bad_shift = self._compute_bad_char_shift(j, mismatch_ch)
                good_shift = self.shift_val[j + 1]
                i += max(bad_shift, good_shift)
        return results

def run_example():
    searcher = BoyerMooreSearcher("ACTG", "ACCA")
    text = "ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC"
    matches = searcher.find_occurrences(text)
    print(matches)

if __name__ == "__main__":
    run_example()
