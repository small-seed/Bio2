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

# Implement the Smith-Waterman dynamic programming for local alignment
def perform_smith_waterman(seq_a, seq_b, sub_matrix, gap_penalty):
    len_a, len_b = len(seq_a), len(seq_b)
    # Create score and traceback matrices
    score_mat = [[0] * (len_b + 1) for _ in range(len_a + 1)]
    trace_mat = [[0] * (len_b + 1) for _ in range(len_a + 1)]
    max_score = 0
    max_row, max_col = 0, 0
    
    # Initialization: first row and column are 0 for local alignment
    # No need to set explicitly as already initialized to 0
    
    # Populate the matrices
    for row in range(1, len_a + 1):
        for col in range(1, len_b + 1):
            diag = score_mat[row-1][col-1] + evaluate_position_score(seq_a[row-1], seq_b[col-1], sub_matrix, gap_penalty)
            vertical = score_mat[row-1][col] + gap_penalty
            horizontal = score_mat[row][col-1] + gap_penalty
            candidates = [diag, vertical, horizontal, 0]
            max_val = max(candidates)
            score_mat[row][col] = max_val
            if max_val == 0:
                trace_mat[row][col] = 0
            else:
                # Determine direction excluding the 0 case
                direction_candidates = [diag, vertical, horizontal]
                trace_mat[row][col] = select_traceback_path(*direction_candidates)
            if max_val > max_score:
                max_score = max_val
                max_row, max_col = row, col
    
    return score_mat, trace_mat, max_score, (max_row, max_col)

# Reconstruct the aligned sequences by following the traceback from the max score position
def reconstruct_local_alignment(trace_matrix, seq_a, seq_b, start_row, start_col):
    aligned_a = []
    aligned_b = []
    row, col = start_row, start_col
    while row > 0 and col > 0 and trace_matrix[row][col] > 0:
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
def run_local_alignment_test():
    sub_matrix = load_substitution_matrix("blosum62.mat")
    sequence_a = "HGWAGWAGG"
    sequence_b = "PHSWGGAGH"
    score_result, trace_result, opt_score, max_pos = perform_smith_waterman(sequence_a, sequence_b, sub_matrix, -8)
    print("Optimal alignment score:", opt_score)
    display_matrix(score_result)
    display_matrix(trace_result)
    alignment_a, alignment_b = reconstruct_local_alignment(trace_result, sequence_a, sequence_b, *max_pos)
    print(alignment_a)
    print(alignment_b)

run_local_alignment_test()
