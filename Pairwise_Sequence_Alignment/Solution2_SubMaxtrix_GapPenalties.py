# Define a function to build a custom scoring matrix for sequence alignments
def build_scoring_matrix(match_score, mismatch_score, symbols):
    scoring_dict = {}
    for symbol1 in symbols:
        for symbol2 in symbols:
            key = symbol1 + symbol2
            scoring_dict[key] = match_score if symbol1 == symbol2 else mismatch_score
    return scoring_dict

# Parse a substitution matrix from an input file
def parse_matrix_file(file_path):
    scoring_dict = {}
    with open(file_path, 'r') as file_handle:
        header_line = file_handle.readline().strip()
        header_tokens = header_line.split('\t')
        num_symbols = len(header_tokens)
        symbols = [token[0] for token in header_tokens]
        
        for row_index in range(num_symbols):
            row_line = file_handle.readline().strip()
            row_tokens = row_line.split('\t')
            for col_index, score_str in enumerate(row_tokens):
                pair_key = symbols[row_index] + symbols[col_index]
                scoring_dict[pair_key] = int(score_str)
    return scoring_dict

# Compute the score for a single aligned position, accounting for gaps
def position_score(char1, char2, scoring_matrix, gap_penalty):
    if char1 == '-' or char2 == '-':
        return gap_penalty
    return scoring_matrix[char1 + char2]

# Compute the total score for a pairwise alignment
def compute_alignment_score(alignment1, alignment2, scoring_matrix, gap_penalty):
    total_score = 0
    for idx in range(len(alignment1)):
        total_score += position_score(alignment1[idx], alignment2[idx], scoring_matrix, gap_penalty)
    return total_score
