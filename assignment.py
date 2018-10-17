import sys
import st

fasta_filename = sys.argv[1]

with open(fasta_filename) as fasta_file:
    for header, seq  in st.get_next_fasta(fasta_file):
        print(header)
        print(len(seq))
