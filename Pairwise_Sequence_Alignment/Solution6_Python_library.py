from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo

def perform_global_alignment(seq1, seq2):
    results = pairwise2.align.globalxx(seq1, seq2)
    return results

def execute_blosum_alignment(seq1, seq2):
    scoring_matrix = MatrixInfo.blosum62
    alignments = pairwise2.align.globalds(seq1, seq2, scoring_matrix, -4, -1)
    for result in alignments:
        print(format_alignment(*result))

def execute_local_alignments(dna_seq1, dna_seq2, prot_seq1, prot_seq2):
    scoring_matrix = MatrixInfo.blosum62
    dna_results = pairwise2.align.localms(dna_seq1, dna_seq2, 3, -2, -3, -3)
    prot_results = pairwise2.align.localds(prot_seq1, prot_seq2, scoring_matrix, -4, -1)
    return dna_results, prot_results

def run_example():
    dna1 = input("Enter first DNA sequence: ")
    dna2 = input("Enter second DNA sequence: ")
    prot1 = input("Enter first protein sequence: ")
    prot2 = input("Enter second protein sequence: ")
    results = execute_local_alignments(dna1, dna2, prot1, prot2)
    print("DNA local alignment:", results[0])
    print("Protein local alignment:", results[1])

if __name__ == "__main__":
    run_example()
