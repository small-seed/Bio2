from solution1_sequence_motif_handler import SequenceMotif
from solution1_sequence_motif_handler import BioSequence
import random

class MotifDiscovery:
    def __init__(self, motif_length=8, sequences=None):
        self.motif_length = motif_length
        self.sequences = sequences if sequences is not None else []
        if self.sequences:
            self.base_alphabet = self.sequences[0].get_alphabet()

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, index):
        return self.sequences[index]

    def get_sequence_length(self, index):
        return len(self.sequences[index])

    def load_sequences_from_file(self, file_path, seq_type):
        with open(file_path, 'r') as file_handle:
            for line in file_handle:
                stripped_line = line.strip().upper()
                if stripped_line:  # Skip empty lines
                    self.sequences.append(BioSequence(stripped_line, seq_type))
        if self.sequences:
            self.base_alphabet = self.sequences[0].get_alphabet()

    def build_motif_from_positions(self, position_list):
        extracted_subseqs = []
        for idx, start_pos in enumerate(position_list):
            subseq_str = self.sequences[idx][start_pos:start_pos + self.motif_length]
            extracted_subseqs.append(BioSequence(subseq_str, self.sequences[idx].sequence_type()))
        return SequenceMotif(extracted_subseqs)

    def compute_additive_score(self, position_list):
        motif_instance = self.build_motif_from_positions(position_list)
        motif_instance._compute_counts()
        count_table = motif_instance.count_matrix
        total_score = 0
        num_positions = len(count_table[0]) if count_table else 0
        num_bases = len(count_table)
        for col in range(num_positions):
            max_in_col = count_table[0][col]
            for row in range(1, num_bases):
                if count_table[row][col] > max_in_col:
                    max_in_col = count_table[row][col]
            total_score += max_in_col
        return total_score

    def compute_product_score(self, position_list):
        motif_instance = self.build_motif_from_positions(position_list)
        motif_instance._compute_weight_matrix()
        prob_table = motif_instance.weight_matrix
        product_score = 1.0
        num_positions = len(prob_table[0]) if prob_table else 0
        num_bases = len(prob_table)
        for col in range(num_positions):
            max_in_col = prob_table[0][col]
            for row in range(1, num_bases):
                if prob_table[row][col] > max_in_col:
                    max_in_col = prob_table[row][col]
            product_score *= max_in_col
        return product_score

    def perform_heuristic_search(self):
        current_positions = [random.randint(0, self.get_sequence_length(i) - self.motif_length)
                             for i in range(len(self.sequences))]
        current_motif = self.build_motif_from_positions(current_positions)
        current_motif._compute_weight_matrix()
        current_score = self.compute_product_score(current_positions)
        best_positions = current_positions[:]
        has_improved = True
        while has_improved:
            for seq_idx in range(len(self.sequences)):
                current_positions[seq_idx] = current_motif.find_best_match_position(self.sequences[seq_idx])
            new_score = self.compute_product_score(current_positions)
            if new_score > current_score:
                current_score = new_score
                best_positions = current_positions[:]
                current_motif = self.build_motif_from_positions(current_positions)
                current_motif._compute_weight_matrix()
            else:
                has_improved = False
        return best_positions

    # Gibbs sampling method
    def execute_gibbs_sampling(self, num_iterations=100):
        initial_positions = [random.randint(0, len(self.sequences[0]) - self.motif_length)
                             for _ in range(len(self.sequences))]
        best_positions = initial_positions[:]
        best_score = self.compute_product_score(initial_positions)
        for _ in range(num_iterations):
            # Select a random sequence to resample
            selected_seq_idx = random.randint(0, len(self.sequences) - 1)
            selected_seq = self.sequences[selected_seq_idx]
            # Temporarily remove the sequence
            temp_positions = initial_positions[:selected_seq_idx] + initial_positions[selected_seq_idx + 1:]
            temp_seqs = self.sequences[:selected_seq_idx] + self.sequences[selected_seq_idx + 1:]
            temp_motif = self.build_motif_from_positions(temp_positions)
            temp_motif._compute_weight_matrix()
            # Compute probabilities for all positions in the selected sequence
            pos_probs = temp_motif.calculate_site_probabilities(selected_seq)
            # Select new position via roulette wheel
            new_pos = self.roulette_wheel_selection(pos_probs)
            initial_positions[selected_seq_idx] = new_pos
            updated_score = self.compute_product_score(initial_positions)
            if updated_score > best_score:
                best_score = updated_score
                best_positions = initial_positions[:]
        return best_positions

    def roulette_wheel_selection(self, probabilities):
        total_sum = sum(0.01 + prob for prob in probabilities)
        random_value = random.random() * total_sum
        cumulative = 0.0
        for idx, prob in enumerate(probabilities):
            cumulative += (prob + 0.01)
            if cumulative >= random_value:
                return idx
        return len(probabilities) - 1
