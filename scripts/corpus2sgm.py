#!/usr/bin/python -tt

'''
File: corpus_to_sgm.py
Date: January 3, 2013
Author: Avneesh Saluja (avneesh@cs.cmu.edu)
Description: simple script that takes a ' ||| ' delimited
file of source-target sentences, and a grammar location,
and writes the location in a per-sentence grammar and writes
out the result to stdout. 
stdin: corpus
arg1: grammar loc
'''

import sys, commands, string

def main():
    grammarLoc = sys.argv[1]
    #include hiero grammar here
    count = 0
    for line in sys.stdin:
        elements = line.strip().split(' ||| ')
        openTag = '<seg grammar="%sgrammar.%d.gz" id="%d">'%(grammarLoc, count, count)
        closeTag = '</seg>'
        divider = ' ||| '
        count += 1
        #print ''.join([openTag, elements[0], closeTag, divider, ' ||| '.join(elements[1:])])
        print ''.join([openTag, elements[0], divider, ' ||| '.join(elements[1:]), closeTag])

if __name__ == "__main__":
    main()
