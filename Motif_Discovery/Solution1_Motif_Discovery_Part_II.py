from Solution1_Motif_Discovery_Part_I import Sequence

class DeterministicMotifFinder:
    """
    A class for deterministic motif discovery in biological sequences.
    Supports exhaustive search, branch-and-bound, and heuristic consensus methods.
    """
    
    def __init__(self, motif_length=8, sequences=None):
        """
        Initialize the motif finder.
        
        :param motif_length: Length of the motif to search for.
        :param sequences: List of Sequence objects (optional).
        """
        self.motif_length = motif_length
        self.sequences = sequences if sequences is not None else []
        self.alphabet = sequences[0].get_alphabet() if sequences else None
    
    def __len__(self):
        """Return the number of sequences."""
        return len(self.sequences)
    
    def __getitem__(self, index):
        """Access a sequence by index."""
        return self.sequences[index]
    
    def sequence_length(self, index):
        """Return the length of the sequence at the given index."""
        return len(self.sequences[index])
    
    def load_sequences_from_file(self, filename, seq_type):
        """
        Load sequences from a file.
        
        :param filename: Path to the file containing sequences (one per line).
        :param seq_type: Type of sequences ("DNA", "RNA", or "PROTEIN").
        """
        self.sequences = []
        with open(filename, 'r') as file:
            for line in file:
                stripped_line = line.strip().upper()
                if stripped_line:  # Skip empty lines
                    self.sequences.append(Sequence(stripped_line, seq_type))
        if self.sequences:
            self.alphabet = self.sequences[0].get_alphabet()
    
    def _create_profile_matrix(self, start_positions):
        """
        Create a profile matrix (count matrix) from start positions in each sequence.
        
        :param start_positions: List of starting indices for the motif in each sequence.
        :return: 2D list where rows are alphabet symbols, columns are motif positions.
        """
        num_sequences = len(self.sequences)
        alphabet_size = len(self.alphabet)
        profile = [[0] * self.motif_length for _ in range(alphabet_size)]
        
        for seq_idx, start_pos in enumerate(start_positions):
            subsequence = self.sequences[seq_idx][start_pos:start_pos + self.motif_length]
            for pos_idx, nucleotide in enumerate(subsequence):
                if nucleotide in self.alphabet:
                    alpha_idx = self.alphabet.index(nucleotide)
                    profile[alpha_idx][pos_idx] += 1
        
        return profile
    
    def calculate_additive_score(self, start_positions):
        """
        Calculate the additive score of a set of start positions (sum of max counts per column).
        
        :param start_positions: List of starting indices.
        :return: The additive score.
        """
        profile = self._create_profile_matrix(start_positions)
        total_score = 0
        for col in range(self.motif_length):
            col_max = max(row[col] for row in profile)
            total_score += col_max
        return total_score
    
    def calculate_multiplicative_score(self, start_positions):
        """
        Calculate the multiplicative score of a set of start positions (product of max counts per column).
        
        :param start_positions: List of starting indices.
        :return: The multiplicative score.
        """
        profile = self._create_profile_matrix(start_positions)
        total_score = 1.0
        for col in range(self.motif_length):
            col_max = max(row[col] for row in profile)
            total_score *= col_max
        return total_score
    
    def _generate_next_position_vector(self, current_positions):
        """
        Generate the next lexicographical position vector.
        
        :param current_positions: Current list of start positions.
        :return: Next position vector or None if at the end.
        """
        next_positions = list(current_positions)
        i = len(next_positions) - 1
        max_positions = [self.sequence_length(idx) - self.motif_length for idx in range(len(next_positions))]
        
        while i >= 0 and next_positions[i] == max_positions[i]:
            i -= 1
        
        if i < 0:
            return None
        
        next_positions[i] += 1
        for j in range(i + 1, len(next_positions)):
            next_positions[j] = 0
        
        return next_positions
    
    def exhaustive_search(self):
        """
        Perform exhaustive search to find the optimal motif positions.
        
        :return: List of optimal start positions.
        """
        if not self.sequences:
            return None
        
        best_score = -1
        best_positions = None
        current_positions = [0] * len(self.sequences)
        
        while current_positions is not None:
            score_val = self.calculate_additive_score(current_positions)
            if score_val > best_score:
                best_score = score_val
                best_positions = list(current_positions)
            current_positions = self._generate_next_position_vector(current_positions)
        
        return best_positions
    
    def _generate_next_vertex(self, current_path):
        """
        Generate the next vertex in the branch-and-bound tree (extend or increment).
        
        :param current_path: Current path in the search tree.
        :return: Next path or None.
        """
        num_sequences = len(self.sequences)
        if len(current_path) < num_sequences:
            # Extend to next sequence
            return current_path + [0]
        else:
            # Increment like next position vector
            return self._generate_next_position_vector(current_path)
    
    def _bypass_subtree(self, current_path):
        """
        Bypass the current subtree by incrementing the last fixed position.
        
        :param current_path: Current partial path.
        :return: Next path at the same level or None.
        """
        # Similar to _generate_next_position_vector but for partial path
        next_path = list(current_path)
        i = len(next_path) - 1
        max_positions = [self.sequence_length(idx) - self.motif_length for idx in range(len(next_path))]
        
        while i >= 0 and next_path[i] == max_positions[i]:
            i -= 1
        
        if i < 0:
            return None
        
        next_path[i] += 1
        for j in range(i + 1, len(next_path)):
            next_path[j] = 0
        
        return next_path
    
    def branch_and_bound_search(self):
        """
        Perform branch-and-bound search to find the optimal motif positions.
        
        :return: List of optimal start positions.
        """
        if not self.sequences:
            return None
        
        best_score = -1
        best_positions = None
        num_sequences = len(self.sequences)
        current_path = []
        
        while current_path is not None:
            if len(current_path) < num_sequences:
                # Internal node: estimate upper bound
                partial_score = self.calculate_additive_score(current_path)
                remaining_contribution = (num_sequences - len(current_path)) * self.motif_length
                upper_bound = partial_score + remaining_contribution
                
                if upper_bound <= best_score:  # Note: <= to handle equality if desired
                    current_path = self._bypass_subtree(current_path)
                else:
                    current_path = self._generate_next_vertex(current_path)
            else:
                # Leaf node: evaluate full solution
                score_val = self.calculate_additive_score(current_path)
                if score_val > best_score:
                    best_score = score_val
                    best_positions = list(current_path)
                current_path = self._generate_next_vertex(current_path)
        
        return best_positions
    
    def consensus_heuristic(self):
        """
        Use the CONSENSUS heuristic to find approximate motif positions.
        Builds progressively by adding one sequence at a time.
        
        :return: List of approximate start positions.
        """
        num_sequences = len(self.sequences)
        if num_sequences < 2:
            return [0] * num_sequences
        
        best_positions = [0] * num_sequences
        max_score = -1
        
        # Step 1: Find best pair for first two sequences
        seq0_max = self.sequence_length(0) - self.motif_length
        seq1_max = self.sequence_length(1) - self.motif_length
        for pos0 in range(seq0_max + 1):
            for pos1 in range(seq1_max + 1):
                partial = [pos0, pos1]
                score_val = self.calculate_additive_score(partial)
                if score_val > max_score:
                    max_score = score_val
                    best_positions[0] = pos0
                    best_positions[1] = pos1
        
        # Step 2: Add remaining sequences one
