from Solution1_Simple_Progressive_MSA_Part_I import PairwiseAligner as PairwiseAlignment
from Solution1_Simple_Progressive_MSA_Part_I import Alignment as MyAlign
from Solution1_Simple_Progressive_MSA_Part_I import Sequence as MySeq
from Solution1_Simple_Progressive_MSA_Part_I import SubstitutionMatrix as SubstMatrix


class ProgressiveMultipleAligner:
    
    def __init__(self, input_sequences, pairwise_aligner: PairwiseAlignment):
        self.input_sequences = input_sequences
        self.aligner = pairwise_aligner
    
    def _align_to_current(self, current_alignment: MyAlign, new_sequence: MySeq):
        num_existing = current_alignment.num_sequences()
        new_aligned_rows = [""] * (num_existing + 1)
        
        # Compute consensus sequence as guide
        consensus_guide = MySeq(current_alignment.consensus(), current_alignment.alignment_type)
        
        # Align consensus to new sequence
        self.aligner.needleman_wunsch(consensus_guide, new_sequence)
        new_pair_align = self.aligner.recover_global_alignment()
        
        source_col = 0
        for pos in range(len(new_pair_align)):
            if new_pair_align[0, pos] == '-':
                # Insert gap in all existing sequences
                for row in range(num_existing):
                    new_aligned_rows[row] += '-'
            else:
                # Copy column from current alignment
                for row in range(num_existing):
                    new_aligned_rows[row] += current_alignment[row, source_col]
                source_col += 1
        
        # Add the aligned new sequence
        new_aligned_rows[num_existing] = new_pair_align[1]
        
        return MyAlign(new_aligned_rows, current_alignment.alignment_type)
    
    def perform_alignment(self):
        if len(self.input_sequences) < 2:
            return MyAlign([str(seq) for seq in self.input_sequences], self.input_sequences[0].get_type() if self.input_sequences else "DNA")
        
        # Initial pairwise alignment of first two sequences
        self.aligner.needleman_wunsch(self.input_sequences[0], self.input_sequences[1])
        current_alignment = self.aligner.recover_global_alignment()
        
        # Iteratively add remaining sequences
        for seq_idx in range(2, len(self.input_sequences)):
            current_alignment = self._align_to_current(current_alignment, self.input_sequences[seq_idx])
        
        return current_alignment


def run_example():
    seq1 = MySeq("ATAGC")
    seq2 = MySeq("AACC")
    seq3 = MySeq("ATGAC")
    
    sub_matrix = SubstMatrix()
    sub_matrix.create_simple_matrix(match_score=1, mismatch_score=-1, alphabet=list("ACGT"))
    
    pair_aligner = PairwiseAlignment(sub_matrix, gap_penalty=-1)
    
    multi_aligner = ProgressiveMultipleAligner([seq1, seq2, seq3], pair_aligner)
    final_alignment = multi_aligner.perform_alignment()
    
    print(final_alignment)
