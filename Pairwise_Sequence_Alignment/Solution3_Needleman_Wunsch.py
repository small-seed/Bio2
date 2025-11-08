# Load a substitution matrix from a tab-delimited file
def load_substitution_matrix(file_path):
    sub_matrix = {}
    with open(file_path, 'r') as input_file:
        header_row = input_file.readline().strip().split('\t')
        num_symbols = len(header_row)
        symbols = [entry[0] for entry in header_row]
        for row_idx in range(num_symbols):
            row_data = input_file.readline().strip().split('\t')
            for col_idx, score_value in enumerate(row_data):
                symbol_pair = symbols[row_idx] + symbols[col_idx]
                sub_matrix[symbol_pair] = int(score_value)
    return sub_matrix

# Evaluate the score for a pair of characters in the alignment, applying gap penalty if needed
def evaluate_position_score(char_a, char_b, sub_matrix, gap_penalty):
    if char_a == '-' or char_b == '-':
        return gap_penalty
    return sub_matrix[char_a + char_b]

# Determine the direction for traceback based on the highest score
def select_traceback_path(diag_val, vertical_val, horizontal_val):
    if diag_val >= vertical_val and diag_val >= horizontal_val:
        return 1  # Diagonal move
    elif vertical_val >= horizontal_val:
        return 2  # Vertical move
    else:
        return 3  # Horizontal move

# Implement the Needleman-Wunsch dynamic programming for global alignment
def perform_needleman_wunsch(seq_a, seq_b, sub_matrix, gap_penalty):
    len_a, len_b = len(seq_a), len(seq_b)
    # Create score and traceback matrices
    score_mat = [[0] * (len_b + 1) for _ in range(len_a + 1)]
    trace_mat = [[0] * (len_b + 1) for _ in range(len_a + 1)]
    
    # Set up the first row (gaps in seq_a)
    for col in range(1, len_b + 1):
        score_mat[0][col] = col * gap_penalty
        trace_mat[0][col] = 3  # Horizontal
    
    # Set up the first column (gaps in seq_b)
    for row in range(1, len_a + 1):
        score_mat[row][0] = row * gap_penalty
        trace_mat[row][0] = 2  # Vertical
    
    # Populate the rest of the matrices
    for row in range(1, len_a + 1):
        for col in range(1, len_b + 1):
            diag = score_mat[row-1][col-1] + evaluate_position_score(seq_a[row-1], seq_b[col-1], sub_matrix, gap_penalty)
            vertical = score_mat[row-1][col] + gap_penalty
            horizontal = score_mat[row][col-1] + gap_penalty
            max_val = max(diag, vertical, horizontal)
            score_mat[row][col] = max_val
            trace_mat[row][col] = select_traceback_path(diag, vertical, horizontal)
    
    return score_mat, trace_mat

# Reconstruct the aligned sequences by following the traceback matrix
def reconstruct_alignment(trace_matrix, seq_a, seq_b):
    aligned_a = []
    aligned_b = []
    row, col = len(seq_a), len(seq_b)
    while row > 0 or col > 0:
        direction = trace_matrix[row][col]
        if direction == 1:
            aligned_a.append(seq_a[row-1])
            aligned_b.append(seq_b[col-1])
            row -= 1
            col -= 1
        elif direction == 2:
            aligned_a.append(seq_a[row-1])
            aligned_b.append('-')
            row -= 1
        else:
            aligned_a.append('-')
            aligned_b.append(seq_b[col-1])
            col -= 1
    return ''.join(reversed(aligned_a)), ''.join(reversed(aligned_b))

# Output a matrix row by row
def display_matrix(mat):
    for row in mat:
        print(row)

# Example execution
def run_global_alignment_test():
    sub_matrix = load_substitution_matrix("blosum62.mat")
    sequence_a = "PHSWG"
    sequence_b = "HGWAG"
    score_result, trace_result = perform_needleman_wunsch(sequence_a, sequence_b, sub_matrix, -8)
    print("Optimal alignment score:", score_result[-1][-1])
    display_matrix(score_result)
    display_matrix(trace_result)
    alignment_a, alignment_b = reconstruct_alignment(trace_result, sequence_a, sequence_b)
    print(alignment_a)
    print(alignment_b)

run_global_alignment_test()
