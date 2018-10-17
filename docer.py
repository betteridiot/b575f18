import sys
from importlib import import_module

def main(file1, file2):
    file1 = import_module(file1)
    module_doc = file1.__doc__
    with open(file2, 'w') as out:
        out.write(module_doc)

if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[2]
    main(infile, outfile)
