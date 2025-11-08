from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def create_multiple_alignment():
    try:
        num_sequences = int(input("Enter the number of protein sequences: "))
        if num_sequences <= 0:
            print("Please enter a positive number of sequences.")
            return
    except ValueError:
        print("Invalid input: Number of sequences must be an integer.")
        return
    
    sequence_records = []
    for index in range(num_sequences):
        seq_input = input(f"Enter sequence {index + 1}: ").strip().upper()
        if not seq_input:
            print(f"Empty sequence provided for {index + 1}. Skipping.")
            continue
        record = SeqRecord(Seq(seq_input), id=f"seq_{index + 1}")
        sequence_records.append(record)
    
    if len(sequence_records) < 2:
        print("Need at least two sequences for alignment.")
        return
    
    alignment_result = MultipleSeqAlignment(sequence_records)
    print("\nGenerated Multiple Sequence Alignment:")
    print(alignment_result)
