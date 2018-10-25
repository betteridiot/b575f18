"""The st2 module provides utility functions for working with nucleic acid and
protein sequence data. This initial form only includes the universal
translation table and a generator for reading fasta sequence files.

Now extended to include more functions.

Attributes:
    transTab1L:     Dictionary where keys are DNA codons and associated
                    values are corresponding amino acid single character
                    codes.
    compStrDNA:     an str translation table that maps complementary bases
                    only for the characters A,C,G,T. Used with str.translate()
"""

#This is the one letter universal translation table
#It handles cases of DNA ambiguity where the encoded amino acid is unambiguous.
#You need to deal with the missing cases where ambiguity codes would result in
#an ambiguous amino acid assignment. It is suggested that you use 'X' in these
#cases as this is the standard character for an unknown amino acid.
#Only Y (pyrimidine), R (purine) and N (any) degeneracy symbols are handled at
#this time. (need to add M,K,W,S,B,D,H,V where appropriate)
#Stop codons are symbolized as X
#Reassign TAA, TAG, TAR and TGA to change the stop codon symbol if desired.

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


def get_next_fasta(fileObject):
    """
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
    """
    
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


def fasta_format(header, seq, linelength = 60):
    """
    usage: fastaFormattedSeq = fasta_format(header, seq, linelength)

    Parameters:
        header: str, the sequence header
        seq: str, the nucleotide or protein sequence, should not contain returns
        linelength: int, default=60, length of the sequence lines
        
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


# =============================================================================
# understand the str.maketrans() function, it can be very useful for
# substituting multiple characters in one function call if the substitution is
# one character to one character. It can also delete a set of characters
# by giving a third argument to maketrans that is a string of the characters
# you would like removed. It is used in conjunction with str.translate().
# Note: you can map multiple characters to the same substitution character:
# maketrans('acgt', 'xxxx') - replace a,c,g or t with x
# =============================================================================
compStrDNA = str.maketrans('ACGTacgt', 'TGCAtgca')
def complement(seq):
    """
    Usage:
        complement(seq)

    Returns the complement of a DNA sequence while preserving the case of the 
    characters. Works only with DNA sequences consisting solely of
    A, C, G, T or N characters. It does not the IUPAC DNA ambiguity codes.

    Args:
        seq (str): a DNA sequence

    Returns:
        (str): the complement of the input DNA sequence with case preserved,
               empty if seq was empty
    """
    return seq.translate(compStrDNA)


def reverse(seq):
    """
    Usage:
        reverse(seq)

    Args:
        seq (str):  intended to be a biological sequence, but may be any string

    Returns:
        (str): the reverse of the input string, empty if seq was empty
    """
    pass


def reverse_complement(seq):
    """
    Usage:
        reverse_complement(seq)

    Returns the reverse complement of the input string representing a DNA 
    sequence. Works only with DNA sequences consisting solely of  A, C, G, T or N 
    characters. Preserves the case of the input sequence.

    Args:
        seq: str, a DNA sequence string

    Returns:
        str,    the reverse complement of the input DNA sequence string. str is
                empty if seq was empty string
    """
    pass


def gen_codon_dict():
    """Generate a dictionary of all codons plus an additional entry ('NNN') for
    codons that contain an ambiguous entry such as N, Y, U etc. All dictionary
    entries are initialized to 0. This dictionary is used to accumulate codon
    counts in the function count_codons.

    Args:
        None
    Returns:
        A dictionary - keys are codons, items are all 0, includes the key "NNN"
        to account for codons that include a symbol not in ['A', 'C', 'G', 'T']

    >>> codonDict = gen_codon_dict()
    >>> sorted(list(codonDict.items()))
    [('AAA', 0), ('AAC', 0), ('AAG', 0), ('AAT', 0), ('ACA', 0), ('ACC', 0), \
    ('ACG', 0), ('ACT', 0), ('AGA', 0), ('AGC', 0), ('AGG', 0), ('AGT', 0), \
    ('ATA', 0), ('ATC', 0), ('ATG', 0), ('ATT', 0), ('CAA', 0), ('CAC', 0), \
    ('CAG', 0), ('CAT', 0), ('CCA', 0), ('CCC', 0), ('CCG', 0), ('CCT', 0), \
    ('CGA', 0), ('CGC', 0), ('CGG', 0), ('CGT', 0), ('CTA', 0), ('CTC', 0), \
    ('CTG', 0), ('CTT', 0), ('GAA', 0), ('GAC', 0), ('GAG', 0), ('GAT', 0), \
    ('GCA', 0), ('GCC', 0), ('GCG', 0), ('GCT', 0), ('GGA', 0), ('GGC', 0), \
    ('GGG', 0), ('GGT', 0), ('GTA', 0), ('GTC', 0), ('GTG', 0), ('GTT', 0), \
    ('NNN', 0), ('TAA', 0), ('TAC', 0), ('TAG', 0), ('TAT', 0), ('TCA', 0), \
    ('TCC', 0), ('TCG', 0), ('TCT', 0), ('TGA', 0), ('TGC', 0), ('TGG', 0), \
    ('TGT', 0), ('TTA', 0), ('TTC', 0), ('TTG', 0), ('TTT', 0)]
    """
    pass


# =============================================================================
# You will need to write the previous function before you can write
# count_codons because you will need a dictionary initialized for all codons
# plus a special "NNN" key to store counts for any codon with ambiguous
# DNA residues
# =============================================================================
def count_codons(seq):
    """
    Usage: 
        count_codons(seq)

    Counts the codons in an input DNA sequence string beginning with the first 
    base, 3 bases at a time (non-overlapping). All codons are represented in the
    dictionary, and if a codon is not seen, its count will be zero. Any codons 
    containing non A,C,G or T characters are counted under the 'NNN' key.

    Args:
        seq (str): DNA sequence string consisting solely of A, C, G, T or N
                   characters, upper or lower case

    Returns:
        (dict): a dictionary containing the counts of all codons in the input
                sequence. Upper case strings representing the codons are the keys,
                and the associated values are the counts. dict is empty if seq was
                empty string
    """
    pass


def count_kmers(seq, kmerLen, resultDict = None):
    """
    Usage:
        count_kmers(seq, kmerLen)

    Counts the occurrence of kmers (substrings) of length kmerLen in a string. 
    kmers overlap by kmerLen-1 residues, that is, each kmer starts one character
    further along the string than the previous kmer. Only complete kmers are counted.

    Args:
        seq (str): input sequence; DNA, RNA or protein
        kmerLen (int): length of the kmers that are counted
        resultDict (dict):  default None, optional dictionary that already contains
                            counts of kmers so that cumulative counts may be obtained
                            for multiple sequences. Common use case: accumulate kmer
                            counts for a sequence and its reverse complement. If None,
                            a new, empty dict is created.

    Returns:
        (dict): dictionary of counts. keys are the sequences of the kmers and
                the values are the integer counts of the number of occurrences
                of those kmers in the sequence. Dictionary is empty if seq was
                empty or shorter than kmerLen
    """
    pass


# =============================================================================
# fill in the doc string for the following function as well as code it.
# it should return the decimal fractional content of Gs + Cs in a DNA sequence.
# Do not consider any characters other than A,C,G,T
# Note: it is only necessary to do this for one strand.
# =============================================================================
def gc_content(seq):
    """

    """
    pass


def translate_dna(seq):
    """
    Usage:
        translate_dna(seq)

    Translates frame 1 of the input DNA sequence, that is, the first codon
    starts at the first base of the input sequence. Uses the universal
    translation code. Can handle IUPAC degeneracy codes that unambiguously
    map to a single amino acid. The translation uses single-character
    amino acid codes. Stop codons are translated as "X". Ambiguous codons
    (not in the translation table) are also represented as X. X is used so
    that the output may be used by sequence search and alignment programs.

    Args:
        seq (str): input sequence; DNA only, upper and/or lower case
    Returns:
        (str): A string of single-letter amino acid codes in upper case. If
                seq is empty, an empty string is returned.
    """
    pass

# =============================================================================
# The next portion is to add a quick summary of the FIRST fasta entry in any
# fasta file. We expect the following (in order):
# 1. Header
# 2. Sequence Length
# 3. The GC content of the first sequence
# 4. The most abundant codon in the first sequence
# 5. The first 10 translated codons

# NOTE: pay VERY CLOSE attention to the expected output shown in the homework page
# =============================================================================


# Add the command line execution heuristics here and conditional imports:
