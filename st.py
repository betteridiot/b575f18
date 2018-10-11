'''The st module provides utility functions for working with nucleic acid and
protein sequence data. This initial form only includes the universal translation
table and a generator for reading fasta sequence files.
'''

#This is the one letter universal translation table
#It handles cases of DNA ambiguity where the encoded amino acid is unambiguous.
#You need to deal with the missing cases where ambiguity codes would result in
#an ambiguous amino acid assignment. It is suggested that you use 'X' in these
#cases as this is the standard character for an unknown amino acid.
#Only Y (pyrimidine), R (purine) and N (any) degeneracy symbols are handled at
#this time. (need to add M,K,W,S,B,D,H,V where appropirate)
#Stop codons are symbolized as X
#Reassign TAA, TAG, TAR and TGA to change the stop codon sybmol if desired.

transTab1L = {
'TTT': 'F', 'TTC': 'F', 'TTY': 'F', 'TTA': 'L', 'TTG': 'L', 'TTR': 'L', 
'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TCN': 'S', 'TCY': 'S', 'TCR': 'S', 
'TAT': 'Y', 'TAC': 'Y', 'TAY': 'Y', 'TAA': 'X', 'TAG': 'X', 'TAR': 'X', 
'TGT': 'C', 'TGC': 'C', 'TGY': 'C', 'TGA': 'X', 'TGG': 'W', 
'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CTY': 'L', 'CTR': 'L', 'CTN': 'L',
    'YTG': 'L', 'YTA': 'L', 
'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCY': 'P', 'CCR': 'P', 'CCN': 'P', 
'CAT': 'H', 'CAC': 'H', 'CAY': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CAR': 'Q', 
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'CGY': 'R', 'CGR': 'R', 'CGN': 'R', 
'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATY': 'I', 'ATG': 'M', 
'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'ACY': 'T', 'ACR': 'T', 'ACN': 'T', 
'AAT': 'N', 'AAC': 'N', 'AAY': 'N', 'AAA': 'K', 'AAG': 'K', 'AAR': 'K', 
'AGT': 'S', 'AGC': 'S', 'AGY': 'S', 'AGA': 'R', 'AGG': 'R', 'AGR': 'R', 
'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GTY': 'V', 'GTR': 'V', 'GTN': 'V', 
'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GCY': 'A', 'GCR': 'A', 'GCN': 'A', 
'GAT': 'D', 'GAC': 'D', 'GAY': 'D', 'GAA': 'E', 'GAG': 'E', 'GAR': 'E', 
'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGY': 'G', 'GGR': 'G', 'GGN': 'G'
}


def get_next_fasta (fileObject):
    '''
usage: for header, seq in get_next_fasta(fileObject):
            process header
            process seq

Parameters:
    fileObject: an open fileObject for a fasta file

Yields:
    header: str, contains the header line of a fasta record
    seq: str, contains the sequence of a fasta record, newlines removed

This is a generator that returns one fasta record's header and
sequence at a time from a multiple fasta file. The return character is
removed from the header. The sequence is returned as one continuous string
with no returns. The returned value is a tuple (header, sequence)
If their is no sequence associated with a header, seq will be an
empty string
Code simplification contributed by Dattatreya Mellacheruvu
01/16/2009, Jeffrey R. de Wet
08/02/2010 refactored to put lines into a temporary list
    '''
    
    header = ''
    seq = ''

    lineList = []
    #The following for loop gets the header of the first fasta
    #record. Skips any leading junk in the file
    for line in fileObject:
        if line.startswith('>'):
            header = line.strip()
            break
    
    for line in fileObject:
        if line.startswith('>'):
            seq = ''.join(lineList)
            lineList = []
            yield header, seq
            header = line.strip()
            seq = ''
        else:
            #seq += line.strip()
            lineList.append(line.strip())
    #yield the last entry
    if header:
        seq = ''.join(lineList)
        yield header, seq


def fasta_format (header, seq, linelength = 60):
    """
usage: fastaFormattedSeq = fasta_format(header, seq, linelength)

Parameters:
    header: str, the sequence header
    seq: str, the nucleotide or protein sequence, should not contain returns
    linelength: int, length of the sequence lines, default=60
    
Returns:
    str, fasta formatted sequence, return added to last line

fasta_format will convert a sequence to a Fasta formated sequence.
The name (definition line) doesn't need to include the > symbol or
return characters. The function will add them. The header and seq (sequence)
should be strings. The sequence should not include return characters

The linelength controls the length of the sequence lines only.
The sequence should not contain return characters.
One of the filters in this module can be used to strip
inappropriate characters before using the formatter.
"""
    
    tempLines = []
    header = header.strip()
    if not header.startswith('>'):
        header = '>' + header
    tempLines.append(header)
    for x in range(0, len(seq), linelength):
        tempLines.append(seq[x:x + linelength])
    #add an empty string so that a \n will be joined at the end
    tempLines.append('')
    return '\n'.join(tempLines)