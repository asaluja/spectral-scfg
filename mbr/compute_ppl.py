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
        logp = float(line.strip().split(':')[1])
        if logp > 0:
            numValid += 1
            corpus_ll += -logp
    ppl = math.exp(-corpus_ll/numValid)
    print "Perplexity: %.3f"%ppl
    print "Reachable sentences: %d/%d"%(numValid, numTotal)
            
            

if __name__ == "__main__":
    main()
