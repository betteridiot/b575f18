#!/usr/bin/env python

"""makeTable.py
    
    Simple program to create a numCols X numRows sized table 
    filled with random float numbers between -1000 and 1000
    
    usage: makeTable.py <numRows> <numCols>
"""
# io is used so we can define the end-of-line character. On Linux/Mac is a \n and windows is \r\n. This is usually unseen by the user.
import csv
import os
import sys
import random


# additional types in the function definition are called typing. Doesn't change anything for running the code, just makes it more transparent to users what that code needs.
def make_row(numCols: int) -> list:
    """Makes a numCols length list of random floats from uniform distribution with parameters a = -1000 and b = 1000
    
    Args:
        numCols (int): number of desired columns for table
    
    Returns:
        row (str): a tab-delimited string representation of row list
    """
    # below is an example of a list comprehension. This automatically creates a list of each element in the for loop all at once.
    row = [random.uniform(-1000, 1000) for col in range(numCols)]
    return row


def make_table() -> int:
    """Driver function for consuming arguments, producing table, and writing output
    
    Returns:
        0 (int): just a signal
    """
    # setting the conditions below allows the program to give feedback as to why the program would fail before it runs
    if len(sys.argv) > 3:
        print("Too many arguments submitted")
    elif len(sys.argv) == 1:
        print("No table dimensions submitted")
    else:
        numRows = int(sys.argv[1])
        numCols = int(sys.argv[2])
        
        # using a context manager (with) that automatically closes the opened file when operation is done
        with open("table.tsv", 'w', newline = '') as outFile:
            tsvWriter = csv.writer(outFile, delimiter = "\t", lineterminator = os.linesep)
            # the underscore is called a throw-away variable, since we don't really care about it or using it.
            for _ in range(numRows):
                tsvWriter.writerow(make_row(numCols))
    return 0

# the section below allows the program to be run from the command line as a program. Furthermore, using this, we can import this script into another program and use the functions we have already typed
if __name__ == "__main__":
    make_table()