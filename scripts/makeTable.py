#!/users/env python
# -*- coding: utf-8 -*-
"""
usage: makeTable.py <numRows> <numCols>
makes a table, table.tsv of dimensions numRows x numCols
containing random numbers, uniform distribution, range -1000 to +1000
"""

import sys
import random

numRows = int(sys.argv[1])
numCols = int(sys.argv[2])

random.seed(0)

outFile = open('table.tsv', 'w')

# Making a function that can make a row at a time
def make_row(numCols):
    '''usage make_rows(<numCols>)
    Returns a list of random float, list length = numCols
    '''
    # the parameters for the uniform distribution
    low = -1000
    high = 1000
    
    # make empty list to capture random numbers
    row = []
    
    # now generate numCols number of random numbers and put them in the list
    for c in range(numCols):
        row.append(random.uniform(low,high))
    return row

# Make a function that makes numRows number of rows
for r in range(numRows):
    row = make_row(numCols)
    
    # Now that we have a row, make each number tab-separated and write them to file
    # the *row is a shortcut that "unpacks" each item of a list
    print(*row, file=outFile, sep='\t')
outFile.close()
