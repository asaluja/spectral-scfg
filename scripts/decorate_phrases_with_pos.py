#!/usr/bin/python -tt
##############################
#File: decorate_phrases_with_pos.py
#Date: September 17, 2013
#Description: script that attaches POS tags to
#phrases, and also modifies the string formats
#to be friendlier with the spectral code. 
#Author: Avneesh Saluja (avneesh@cs.cmu.edu)
##############################

import sys, commands, string

def main():
    posFile = open(sys.argv[1], 'r')
    posLines = posFile.read().splitlines()
    counter = 0
    for line in sys.stdin:
        phrases = line.strip().split(' | ')
        phrasesPOS = posLines[counter].strip().split(' | ')
        if len(phrases) != len(phrasesPOS):
            sys.stderr.write("Problem: number of phrases in segmented phrases file and segmented POS file do not match.  Skipping sentence.\n");
        else:
            line2print = []
            for i, phrase in enumerate(phrases):
                phrase_reformat = '_'.join(phrase.split()) #replace whitespace with '_'
                POS_reformat = '_'.join(phrasesPOS[i].split())
                line2print.append('#%'.join([phrase_reformat, POS_reformat]))
            print '|'.join(line2print)
        counter += 1

if __name__ == "__main__":
    main()
