import sys
import os

# Top-level commands
print('Top-level loaded')
print(sys.argv[1])

# Encapsulating operations within functions
def foo():
    print('Marcus was here')
    return None

# What happens from the command line
def main(bar):
    print('Command line curator')
    print(bar)
    return None

if __name__ == '__main__':
    print('CLI')
    main(sys.argv[1])

