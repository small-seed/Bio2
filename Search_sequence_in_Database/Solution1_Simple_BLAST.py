# Function to load sequences from a file into a list
def load_sequences_from_file(file_path):
    sequences = []
    with open(file_path, 'r') as file_handle:
        for line in file_handle:
            sequences.append(line.strip())
    return sequences

# Function to create a hash table of k-mers from the query sequence
def create_kmer_index(query_sequence, kmer_length):
    kmer_to_positions = {}
    query_len = len(query_sequence)
    for start_idx in range(query_len - kmer_length + 1):
        kmer = query_sequence[start_idx:start_idx + kmer_length]
        if kmer in kmer_to_positions:
            kmer_to_positions[kmer].append(start_idx)
        else:
            kmer_to_positions[kmer] = [start_idx]
    return kmer_to_positions

# Function to identify exact k-mer matches between a sequence and the query index
def find_kmer_hits(target_sequence, kmer_index, kmer_length):
    matches = []
    target_len = len(target_sequence)
    for start_idx in range(target_len - kmer_length + 1):
        kmer = target_sequence[start_idx:start_idx + kmer_length]
        if kmer in kmer_index:
            for query_start in kmer_index[kmer]:
                matches.append((query_start, start_idx))
    return matches

# Function to extend a seed match bidirectionally using a 50% identity threshold
def extend_seed_match(target_sequence, seed_hit, query_sequence, kmer_length):
    query_start, target_start = seed_hit
    forward_matches = 0
    forward_length = 0
    max_forward = 0
    while (query_start + kmer_length + forward_length < len(query_sequence) and
           target_start + kmer_length + forward_length < len(target_sequence)):
        if query_sequence[query_start + kmer_length + forward_length] == \
           target_sequence[target_start + kmer_length + forward_length]:
            forward_matches += 1
        forward_length += 1
        if 2 * forward_matches >= forward_length:
            max_forward = forward_length
        else:
            break
    
    backward_matches = 0
    backward_length = 0
    max_backward = 0
    while (query_start > backward_length and target_start > backward_length):
        if query_sequence[query_start - 1 - backward_length] == \
           target_sequence[target_start - 1 - backward_length]:
            backward_matches += 1
        backward_length += 1
        if 2 * backward_matches >= backward_length:
            max_backward = backward_length
        else:
            break
    total_length = kmer_length + max_forward + max_backward
    total_matches = kmer_length + forward_matches + backward_matches  # Assuming seed is perfect match
    adjusted_query_start = query_start - max_backward
    adjusted_target_start = target_start - max_backward
    return (adjusted_query_start, adjusted_target_start, total_length, total_matches)

# Function to find the best extended hit for a single target sequence
def find_best_hit_for_sequence(target_sequence, query_sequence, kmer_index, kmer_length):
    seed_hits = find_kmer_hits(target_sequence, kmer_index, kmer_length)
    if not seed_hits:
        return ()
    best_score = -1
    best_alignment = ()
    for hit in seed_hits:
        extended = extend_seed_match(target_sequence, hit, query_sequence, kmer_length)
        current_score = extended[3]
        if (current_score > best_score or
            (current_score == best_score and extended[2] < best_alignment[2])):
            best_score = current_score
            best_alignment = extended
    return best_alignment

# Main function to search a database for the best alignment to the query
def search_database_for_best_alignment(sequence_database, query_sequence, kmer_length):
    kmer_index = create_kmer_index(query_sequence, kmer_length)
    global_best_score = -1
    global_best = (0, 0, 0, 0, 0)
    for db_index, db_sequence in enumerate(sequence_database):
        candidate = find_best_hit_for_sequence(db_sequence, query_sequence, kmer_index, kmer_length)
        if candidate:
            candidate_score = candidate[3]
            if (candidate_score > global_best_score or
                (candidate_score == global_best_score and candidate[2] < global_best[2])):
                global_best_score = candidate_score
                global_best = (*candidate, db_index)
    return global_best if global_best_score >= 0 else ()
