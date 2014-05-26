#!/usr/bin/python -tt

import sys, commands, string, math

'''
File: compute_ppl.py
Date: May 16, 2014
Description: given a file that is the output of force-decode-all.sh,
this script goes through the file and computes per sentence perplexity.
The final result is the per sentence perplexity of the reference translation
for a certain number of sentences in the training set. 
'''

def main():
    corpus_ll = 0
    numValid = 0
    numTotal = 0
    for line in sys.stdin:
        numTotal += 1
        vals = line.strip().split(': ')
        if vals[1] != '':
            logp = float(vals[1].split('|')[0])
            normalizer = float(vals[1].split('|')[1])
            normalized = logp - normalizer
            if logp > 0:
                numValid += 1
                corpus_ll += normalized
    ppl = math.exp(-corpus_ll/numValid)
    print "Conditional LL of reference: %.3f"%corpus_ll
    print "Perplexity: %.5g"%ppl
    print "Reachable sentences: %d/%d"%(numValid, numTotal)


if __name__ == "__main__":
    main()
