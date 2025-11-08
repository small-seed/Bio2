from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import re

def perform_nucleotide_blast(fasta_file_path):
    sequence_record = SeqIO.read(open(fasta_file_path), "fasta")
    
    blast_handle = NCBIWWW.qblast("blastn", "nt", sequence_record.format("fasta"))
    
    with open("blast_output.xml", "w") as output_file:
        output_file.write(blast_handle.read())
    
    blast_handle.close()

def analyze_blast_results(xml_file_path, e_value_cutoff=0.001):
    with open(xml_file_path) as xml_handle:
        blast_records = NCBIXML.parse(xml_handle)

    blast_record = next(blast_records)
    
    print("*** Significant Alignments (E-value < {:.3f}) ***".format(e_value_cutoff))
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < e_value_cutoff:
                print("**** Alignment ****")
                print("Subject sequence:", alignment.title)
                print("Alignment length:", alignment.length)
                print("E-value:", hsp.expect)
                print(hsp.query[:75] + "...")
                print(hsp.match[:75] + "...")
                print(hsp.sbjct[:75] + "...")
    
    if blast_record.alignments:
        first_align = blast_record.alignments[0]
        print("\nFIRST ALIGNMENT DETAILS:")
        print("Accession ID:", first_align.accession)
        print("Hit ID:", first_align.hit_id)
        print("Description:", first_align.hit_def)
        print("Alignment length:", first_align.length)
        print("Number of HSPs:", len(first_align.hsps))
        
        first_hsp = first_align.hsps[0]
        print("Best HSP E-value:", first_hsp.expect)
        print("Best HSP Score:", first_hsp.score)
        print("Best HSP Length:", first_hsp.align_length)
        print("Best HSP Identities:", first_hsp.identities)
        print("Best HSP Alignment:")
        print("Query:", first_hsp.query)
        print("Match:", first_hsp.match)
        print("Subject:", first_hsp.sbjct)
    
    print("\nTOP 10 ALIGNMENTS:")
    for idx in range(min(10, len(blast_record.alignments))):
        align = blast_record.alignments[idx]
        print(f"\n--- Alignment {idx + 1} ---")
        print("Accession ID:", align.accession)
        print("Description:", align.hit_def)
        for hsp in align.hsps:
            print("  HSP E-value:", hsp.expect)
    
    print("\nORGANISMS IN TOP 20 ALIGNMENTS:")
    organism_list = []
    for idx in range(min(20, len(blast_record.alignments))):
        align = blast_record.alignments[idx]
        description = align.hit_def
        match = re.search(r"\[([^\]]+)\]", description)
        if match:
            organism_list.append(match.group(1))
    
    for organism in organism_list:
        print(organism)
