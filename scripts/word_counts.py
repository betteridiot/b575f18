#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 11:20:27 2018

@author: jrdewet
"""

import string
import sys



def count_words(wordCounts, words):
    '''
    '''
    for word in words:
        wordCounts[word] = wordCounts.get(word, 0) + 1
    return wordCounts

def remove_punc(line):
    '''
    '''
    for punc in string.punctuation:
        line  = line.replace(punc, ' ')
    return line

def sort_word_counts(wordCounts):
    '''
    Arguments:
        wordCounts - dictionary, words are keys, values are counts
    Returns:
        list of key-value pairs sorted so that the words are grouped
        (sorted) by count and are alphabetized within each count group.
    '''
    return sorted(wordCounts.items(), key=lambda item: (item[1], item[0]))
    
    
def main():
    wordCounts = dict()
    
    fileName = sys.argv[1]
    infile = open(fileName, encoding='utf8')
    
    #get to the actual start of the book
    for line in infile:
        if line.startswith('*** START'):
            break
    
    #process the remaining text
    for line in infile:
        #Skip the chapter and book lines
        if line.startswith('CHAPTER') or line.startswith('BOOK'):
            continue
        line = line.lower()
        line = remove_punc(line)
        words = line.split()
        wordCounts = count_words(wordCounts, words)
    infile.close()
    wordCounts = list(sort_word_counts(wordCounts))
    
    print(wordCounts[:50])
    print(wordCounts.reverse()[:50])
    ###############################################################
    ########## You need to get the sorting done after this  #######
    ########## and then print the 50 least frequent and 50  #######
    ########## most frequent words and their associated     #######
    ########## counts. Print each of these as a list of     #######
    ########## 2-tuples (word, count)                       #######
    ###############################################################
    
    
if __name__ == '__main__':
    main()
