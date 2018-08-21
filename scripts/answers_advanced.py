#!/usr/bin/env python
import os
import sys
import typing
import traceback


# Simple type alias that means a value can be integer or float
Number = typing.TypeVar('Number', int, float)


def row_mean(row: typing.List[Number]) -> float:
    """ Create a function that gets the mean of the row
    
    Args:
        row (typing.List[Number]): a list of numbers, may be integers or floats
        
    Returns:
        float: the mean of the row
    """
    return sum(row)/len(row)


def online_column_means(col_entries: typing.List[Number], current_means: float = 0, n: int = 0) -> float:
    """ Computes the online means for all given columns.

    An online mean is a formula that updates a mean when presented with 
    an additional value. This allows one to compute the mean as values are
    presented instead of loading all of the values at the same time.
    
    Args:
        col_entries (Number): an integer or float used to reevalute current_mean
        current_mean (float): current column mean
        n (int): number of column entries not-including new value
        
    Returns:
        float: new column mean after appending new entry
    """
    if current_means is None:
        current_means = [0] * len(col_entries)
    n += 1
    for c in range(len(col_entries)):
        delta = col_entries[c] - current_means[c]
        current_means[c] = current_means[c] + delta / n
    return (current_means, n)


def main(path_to_array: str) -> int:
    """ The main driver function of the script. All top-level code is put here.
    
    When the `if __name__ == '__main__'` condition is used, always call your 
    main function `main()` as per convention. Any other function can be named
    whatever you want.
    
    Args:
        array (str): a file path to a tab-delimited (n x m) rectangular array file
    
    Returns:
        int: signal status. 0 if success 1 if failure
    
    Raises:
        FileNotFoundError: invalid path to array file
    """
    if os.path.isfile(path_to_array):
        try:
            with open(path_to_array, 'r') as array:
                count = 0
                col_means = None
                for row in array:
                    columns = [float(col) for col in row.split()]
                    (col_means, count) = online_column_means(columns, col_means, count)
                    print(col_means, row_mean(columns))
            return 0
        except Exception:
            traceback.print_exception(*sys.exc_info())
            return -1
    else:
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