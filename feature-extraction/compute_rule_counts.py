#!/usr/bin/python -tt

import sys, commands, string

def addToDict(newRule, dictCounts):
    if newRule not in dictCounts:
        dictCounts[newRule] = 1
    else:
        dictCounts[newRule] += 1

def printToFile(counts, filename):
    fh = open(filename, 'w')
    print >> fh, "Number of Rules: %d"%(len(counts))
    d_view = [ (v, k) for k, v in counts.iteritems() ]
    d_view.sort(reverse=True)
    sumCounts = 0
    for v,k in d_view:
        print >> fh, "%s: %d"%(k, v)
        sumCounts += v
    print >> fh, "Count of rules: %d"%sumCounts
    fh.close()

inTdictCounts = {}
preTdictCounts = {}
for line in sys.stdin:
    ruleID = int(line.strip().split(' [')[0])
    tgtPhraseAndInfo = line.strip().split(' ||| ')[2]
    rule = line.strip().split(' ||| ')
    NT = rule[0].split()[1]
    if "In:" in tgtPhraseAndInfo: #in-term
        tgtSide = rule[2][:rule[2].index('In:')].rstrip()
        newRule = ' ||| '.join([NT, rule[1], tgtSide])
        addToDict(newRule, inTdictCounts)
    else:
        if "Out:" in tgtPhraseAndInfo:
            tgtSide = rule[2][:rule[2].index('Out:')].rstrip()
            newRule = ' ||| '.join([NT, rule[1], tgtSide])
            addToDict(newRule, preTdictCounts)
printToFile(inTdictCounts, sys.argv[1])
printToFile(preTdictCounts, sys.argv[2])

