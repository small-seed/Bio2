from Bio.Seq import Seq
from Bio import motifs

# Define a list of sequence instances for the motif
seq_list = [
    Seq("TATAA"),
    Seq("TATTA"),
    Seq("TTTAT"),
    Seq("TATAC")
]

# Create the motif object from the sequences
motif_obj = motifs.create(seq_list)

# Output basic information about the motif
print(type(motif_obj))
print(motif_obj)
print(len(motif_obj))
print(motif_obj.consensus)
print(motif_obj.pwm)
print(motif_obj.counts)
print(motif_obj.pssm)

# Generate a weblogo image for visualization
motif_obj.weblogo("custom_motif_logo.png")

# Normalize the counts matrix with pseudocounts
normalized_pwm = motif_obj.counts.normalize(pseudocounts=0.5)

# Compute the log-odds PSSM from the normalized PWM
log_odds_pssm = normalized_pwm.log_odds()

# Display the normalized PWM and PSSM
print(normalized_pwm)
print(log_odds_pssm)

# Define a test sequence to search within
target_sequence = Seq("TTTTATACACTGCATATAACAACCCAAGCATTATAA")

# Find exact matches of the motif instances in the test sequence
for start_pos, matched_seq in motif_obj.instances.search(target_sequence):
    print(start_pos, " ", matched_seq)

# Search for matches using the PSSM with a score threshold
for start_position, match_score in log_odds_pssm.search(target_sequence, threshold=4.0):
    print("Position %d: score = %5.3f" % (start_position, match_score))

# Calculate PSSM scores for every possible position in the sequence
print(log_odds_pssm.calculate(target_sequence))
