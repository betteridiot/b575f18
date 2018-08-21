#!/usr/bin/env python
""" Usage: python <uniqname>_average.py <path/to/array.dat>

The purpose of this program is to process an array-like data structure and
iteratively determine an online/running mean of each column, as well as 
report the current row's mean.
"""
import os
import sys


def row_mean(row):
    """ Create a function that gets the mean of the row
    
    Try starting on this function first. Once you get the appropriate 
    output here, move on to online_column_means function.
    
    Args:
        row (list): a list of numbers, may be integers or floats
        
    Returns:
        float: the mean of the row
    """
    return round(sum(row)/len(row), 4)


def online_column_means(col_entries, n, current_means = None):
    """ Computes the online means for all given columns.

    An online mean is a formula that updates a mean when presented with 
    an additional value. This allows one to compute the mean as values are
    presented instead of loading all of the values at the same time.
    
    Args:
        col_entries (list): an integer or float used to reevaluate current_means
        n (int): nth row
        current_means (list): current column means

    Returns:
        list: new column means after appending new entries
    """
    if current_means is None:
        current_means = [0] * len(col_entries)
    for i, col in enumerate(col_entries):
        delta = col - current_means[i]
        current_means[i] = round(current_means[i] + (delta / n), 4)
    return current_means


def main(path_to_array):
    """ The main driver function of the script. All top-level code is put here.
    
    When the `if __name__ == '__main__'` condition is used, always call your 
    main function `main()` as per convention. Any other function can be named
    whatever you want.
    
    Args:
        path_to_array (str): a file path to a tab-delimited (n x m) rectangular array file
    
    Returns:
        int: signal status. 0 if success 1 if failure
    
    Raises:
        FileNotFoundError: invalid path to array file
    """
    if os.path.isfile(path_to_array):
        """Below, you will find the `with` statement. This is called a context
        manager. They are often used to automatically open and close your
        files when you are done working with them.
        
        with open('filename', 'r') as infile:
            do_something(infile)
        
        is the equivalent to:
        
        infile = open('filename','r')
        do_something(infile)
        infile.close()
        """
        with open(path_to_array, 'r') as array:
            count = 0
            col_means = None
            for row in array:
                count += 1
                columns = [float(col) for col in row.strip().split("\t")]
                col_means = online_column_means(columns, count, col_means)
                print(col_means, row_mean(columns))
        return 0
    else:
        # Raise is a way for users to define custom errors that
        # are less ambiguous than the built-in errors.
        raise FileNotFoundError("Invalid path to array file")
        return -1


""" The conditional below is used for writing a script that expects to be
re-used. The utility of this allows a script to have 2 main uses: 1) if the 
script is run from the command line as `python my_script.py`, it will run as
expected. 2) The script can now be imported into another script without running
when it is imported. '__main__' is the name of the scope in which top-level code
executes. A module can discover whether or not it is running in the main scope 
by checking its own __name__, which allows a common idiom for conditionally 
executing code in a module when it is run as a script or with python -m but not 
when it is imported
"""
if __name__ == "__main__":
    try:
        if main(sys.argv[1]) == 0:
            print('Complete')
        else:
            print('Failure')
    except IndexError:
        print('No path to array file submitted')