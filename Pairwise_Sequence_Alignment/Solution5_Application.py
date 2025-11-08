def compute_longest_common_subsequence(seq_a, seq_b, symbols="ACGT"):
    scoring_matrix = build_scoring_matrix(1, 0, symbols)
    score_matrix, trace_matrix = perform_needleman_wunsch(seq_a, seq_b, scoring_matrix, 0)
    aligned_a, aligned_b = reconstruct_alignment(trace_matrix, seq_a, seq_b)
    common_subseq = ''.join([aligned_a[i] for i in range(len(aligned_a)) if aligned_a[i] == aligned_b[i]])
    return common_subseq

def calculate_edit_distance(seq_a, seq_b, symbols="ACGT"):
    scoring_matrix = build_scoring_matrix(0, -1, symbols)
    score_matrix, trace_matrix = perform_needleman_wunsch(seq_a, seq_b, scoring_matrix, -1)
    edit_dist = -score_matrix[-1][-1]
    return edit_dist

def find_longest_common_substring(seq_a, seq_b, symbols="ACGT"):
    max_length = max(len(seq_a), len(seq_b))
    large_penalty = -(max_length + 1)
    scoring_matrix = build_scoring_matrix(1, large_penalty, symbols)
    score_matrix, trace_matrix, _, max_pos = perform_smith_waterman(seq_a, seq_b, scoring_matrix, large_penalty)
    aligned_a, aligned_b = reconstruct_local_alignment(trace_matrix, seq_a, seq_b, *max_pos)
    return aligned_a
