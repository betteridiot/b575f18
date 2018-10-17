"""Short description of program

Longer description of program

@author: mdsherm@umich.edu
"""

import sys

def pass_class(student):
    return f'{student} passes with an A++'


def main(filename):
    """
    """
    with open(filename) as infile:
        for line in infile:
            print(line.strip())
    return None

if __name__ == '__main__':
    file = sys.argv[1]
    main(file)
