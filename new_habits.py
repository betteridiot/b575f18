def foo():
    print(__name__)
    pass

# Only happens at the command line
if __name__ == '__main__':
    # At command line is the only time I
    # care about command line args
    import sys
    filename = sys.argv[1]
    
    with open(filename) as ff:
        header, sequence = next(get_next_fasta(ff))
        print('First header:', header)
        print('First seqlen:', len(sequence))
        print('GC..', gc_content(sequence))
        
        print('Trna....', translate_dna(sequence)[:10])