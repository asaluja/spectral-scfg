#!/usr/bin/python -tt

'''
File: filter_unreachable_training.py
Description: this script takes the STDERR of the tree_to_rule.py script 
as an argument, and the parallel training corpus as STDIN and removes
the sentences that are unreachable by the grammar.  Sentences are unreachable
if the minimal grammar extracts a rule with arity > 2 from a sentence. 
arg1: STDERR of tree_to_rule.py script
STDIN: parallel sentence corpus
'''

import sys, commands, string

def processInvalidSentences(filename):
    invalidSentences = {}
    for line in open(filename, 'r'):
        sentence_id = line.strip().split(':')[0]
        sentNum = int(sentence_id.split()[1])
        invalidSentences[sentNum] = 1
    return invalidSentences

def main():
    invalidSentences = processInvalidSentences(sys.argv[1])
    for lineNum,line in enumerate(sys.stdin):
        if lineNum not in invalidSentences:
            print line.strip()

if __name__ == "__main__":
    main()

