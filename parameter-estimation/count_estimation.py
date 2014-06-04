#!/usr/bin/python -tt

'''
File: count_estimation.py
Date: May 25, 2014
Description: simple file that, given a 'full' directory (namely, 
a derivation tree for every sentence pair with the rules for every
NT written out), estimates the counts for every rule.  In this   code, 
we filter all rules with number of NTs > 2 on the RHS.  
Update: 28 May: Removed addone functionality, shifted it to featurize_rules.py
'''

import sys, commands, string, cPickle, os, gzip, re, getopt
from tree import tree

Counts = {}

def CheckArity(rule):
    srcRule = rule.split(' ||| ')[1]
    return len(re.findall(r'\[([^]]+)\]', srcRule))

def ProcessInvalidRules(filename):
    invalid_rules = {}
    for line in open(filename, 'rb'):
        sentNum_Str = line.strip().split(':')[0]
        sentNum = int(sentNum_Str.split()[1])
        invalid_rules[sentNum] = 1
    return invalid_rules

'''
This function handles rule formats of the type '[0] ||| [1,1] der [2,2] ||| [2,2] the [1,1]' etc.
'''
def IncrementCounts(training_tree):
    rule = training_tree.rule
    arity = CheckArity(rule)
    if arity < 3:
        Counts[rule] = Counts[rule] + 1 if rule in Counts else 1
    else:
        sys.stderr.write("Rule %s has more than 2 NTs, not including in counts\n"%rule)
    for child in training_tree.children:
        IncrementCounts(child)

def ReformatRule(rule):
    elements = rule.split(' ||| ')
    expr = re.compile(r'\[([^]]+)\]')
    new_src = []
    NT_counter = 1
    for item in elements[1].split():
        if expr.match(item):
            new_src.append("[X,%d]"%NT_counter)
            NT_counter += 1
        else:
            new_src.append(item)
    new_src_str = ' '.join(new_src)
    new_tgt = []
    for item in elements[2].split():
        if expr.match(item):
            NT_counter = int(re.findall(expr, item)[0])
            new_tgt.append("[X,%d]"%NT_counter)
        else:
            new_tgt.append(item)
    new_tgt_str = ' '.join(new_tgt)
    return ' ||| '.join([elements[0], new_src_str, new_tgt_str])
        
'''
this function handles rules of the type '[X] ||| [X] [X] ||| [1] [2]', etc.
'''
def IncrementCountsByLine(filehandle):
    for rule in filehandle:
        arity = CheckArity(rule.strip())
        if arity < 3:
            rule_formatted = ReformatRule(rule.strip())
            Counts[rule_formatted] = Counts[rule_formatted] + 1 if rule in Counts else 1
        else:
            sys.stderr.write("Rule %s (unformatted) has more than 2 NTs, not including in counts\n"%rule)

def WriteOutCounts(fileout):
    reformatted = {}
    binary = 0
    unary = 0
    preterm = 0
    for rule in Counts:
        arity = CheckArity(rule) #first, keep track of the counts of rules across pre-term, unary, and binary
        if arity == 0: 
            preterm += 1
        elif arity == 1:
            unary += 1
        elif arity == 2:
            binary += 1
        else:
            sys.stderr.write("Error! Writing out rule with # of NTs on RHS > 2!\n")
        srcKey = ' ||| '.join(rule.split(' ||| ')[:-1])
        tgtKey = rule.split(' ||| ')[-1]
        srcDict = reformatted[srcKey] if srcKey in reformatted else {}
        srcDict[tgtKey] = Counts[rule]
        reformatted[srcKey] = srcDict
    cPickle.dump(reformatted, open(fileout, "wb"))
    print "Number of pre-terminal rules (including singletons): %d"%preterm
    print "Number of unary rules: %d"%unary
    print "Number of binary rules: %d"%binary

def main():
    (opts, args) = getopt.getopt(sys.argv[1:], 'n:')
    invalid_rules = None
    for opt in opts:
        if opt[0] == '-n':
            invalid_rules = ProcessInvalidRules(opt[1])
    minrule_grammars_loc = args[0]
    fileout = args[1]
    numSentences = len(os.listdir(minrule_grammars_loc))
    for lineNum in range(0, numSentences):
        minrule_fh = gzip.open(minrule_grammars_loc + "grammar.%d.gz"%lineNum)
        if invalid_rules is not None and lineNum in invalid_rules:
            IncrementCountsByLine(minrule_fh)
        else: #either invalid rules is none (in which case we don't need it) or rule is not invalid
            sync_tree = tree(0, None, None, minrule_fh)
            IncrementCounts(sync_tree)
    WriteOutCounts(fileout)
    print "Counts complete"

if __name__ == "__main__":
    main()
