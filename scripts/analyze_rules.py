#!/usr/bin/python -tt

'''
File: analyze_rules.py
Date: December 24, 2013
Author: Avneesh Saluja (avneesh@cs.cmu.edu)
Description: given a comprehensive list of rules, provides some summary statistics.
'''

import sys, commands, string, re

preTerms = {}
inTerms = {}
invalidRules = {}

def checkArity(rule):
    return len(re.findall(r'\[([^]]+)\]', rule))

for line in sys.stdin:
    arity = checkArity(line.strip().split(' ||| ')[1])
    if arity == 0:
        if line.strip() not in preTerms: 
            preTerms[line.strip()] = 1
    elif arity == 1 or arity == 2:
        if line.strip() not in inTerms:
            inTerms[line.strip()] = 1
    else:
        if line.strip() not in invalidRules:
            invalidRules[line.strip()] = 1

print "Total number of unique rules: %d"%(len(preTerms) + len(inTerms) + len(invalidRules))
print "Total number of pre-terminal rules: %d"%(len(preTerms))
print "Total number of in-terminal rules: %d"%(len(inTerms) + len(invalidRules))
print "Out of in-terminal rules, total number that are not valid: %d"%(len(invalidRules))
