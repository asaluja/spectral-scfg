#!/usr/bin/python -tt

'''
File: tree_to_rule.py
Date: September 30, 2013
Description: this script takes the output of the "RC"
code (an implementation of a linear-time shift-reduce based
algorithm to extract "minimal" synchronous grammar rules
from word alignments.  See the following paper for more details:
H. Zhang, D. Gildea, and D. Chiang: Extracting Synchronous Grammar Rules From Word-Level Alignments
in Linear Time (COLING 2008)) and converts it to a list of the rules used. 
arg1: source-target sentence pairs, one line per pair, each pair separated by ' ||| '
arg2: alignments for each sentence pair
stdin: output of rc code
stdout: list of rules
option/flag 1: -d --> debug, prints out full output (if flag off, prints out Hiero output)
option/flag 2: -z --> print out to .gz file: need to provide directory where grammars are written, e.g.:
-z/usr0/home/avneesh/spectral-scfg/data
#Author: Avneesh Saluja (avneesh@cs.cmu.edu)
'''

import os, sys, commands, string, collections, re, getopt, gzip
NTLabel = 0

def createAlignmentDicts(lines):
    linesDicts = []
    for line in lines:
        lineDict = collections.defaultdict(list)
        for alignment in line.strip().split():
            k,v = tuple([int(idx) for idx in alignment.split("-")])
            lineDict[k].append(v)
        linesDicts.append(lineDict)
    return linesDicts

def findEndIdx(substring):
    counter = 0
    idx = 0
    while (True):
        char = substring[idx]
        if char == '(':
            counter += 1
        elif char == ')':
            counter -= 1
        if counter == 0:
             break
        idx += 1
    return idx

def extractTgtYieldOfChild(substring, alignDict):
    srcIdxs = [idx - 1 for idx in map(int, substring.replace('(', ' ').replace(')', ' ').split())]
    tgtIdxs = []
    for srcIdx in srcIdxs:
        tgtIdxs.extend(alignDict[srcIdx]) 
    minTgt = min(tgtIdxs)
    return minTgt
            
def extractRules(count, original, tree, src, tgt, alignDict, nodeLabel, zipOut, sentFile = None):
    global NTLabel
    ruleSrcStr = []
    ruleTgtStr = []
    children = []
    idx = 0
    srcNTCounter = 0
    isLex = False
    while (idx < len(tree)):
        if tree[idx] == '(' and idx == 0: #opening (: just skip
            idx += 1
        elif tree[idx] == ')' and idx == len(tree) - 1: #last ')' in expression: just skip; should break after this            
            idx += 1
        elif tree[idx] == " ": #ignore white space
            idx += 1
        elif tree[idx] == '(' and idx > 0: #new child
            srcNTCounter += 1
            NTLabel += 1 #encountered a new child, so increment NTLabel
            endChildIdx = findEndIdx(tree[idx:]) + idx            
            tgtYieldIdx = extractTgtYieldOfChild(tree[idx:(endChildIdx+1)], alignDict)
            children.append((tree[idx:(endChildIdx+1)], NTLabel)) #append the s-expr corresponding to child, and also child NT label
            recentNT = NTLabel if original else 'X'
            ruleSrcStr.append("[%s,%d]"%(recentNT, srcNTCounter)) #append NT symbol plus corresponding src NT index
            ruleTgtStr.append(("[%s,%d]"%(recentNT, srcNTCounter), tgtYieldIdx)) #append same info to tgt, but in addition the tgt idx for downstream sorting
            idx = endChildIdx + 1 #update to end of child
        else: #we have a number: replace the number with the appropriate string
            isLex = True
            strNum = ""
            while (tree[idx] != " "): #read in entire number
                strNum += tree[idx]
                idx += 1
            srcIdx = int(strNum) - 1
            ruleSrcStr.append(src[srcIdx]) #append token
            tgtIdxs = alignDict[srcIdx]
            for tgtIdx in tgtIdxs: #due to multiple alignments, can append multiple tokens
                ruleTgtStr.append((tgt[tgtIdx], tgtIdx))
    ruleTgtStr = list(set(ruleTgtStr)) #make them unique
    ruleTgtStrSorted = [tup[0] for tup in sorted(ruleTgtStr, key = lambda tup: tup[1])] #sort by target order
    srcStr = ' '.join(ruleSrcStr)
    tgtStr = ' '.join(ruleTgtStrSorted)
    phrase_pair = "%s ||| %s"%(srcStr, tgtStr)
    printRule(count, srcNTCounter, original, isLex, nodeLabel, phrase_pair, zipOut, sentFile)
    for child in children: #recursively do the same for children
        extractRules(count, original, child[0], src, tgt, alignDict, child[1], zipOut, sentFile)

def printRule(count, srcNTCounter, original, isLex, nodeLabel, phrase_pair, zipOut, sentFile=None):
    if zipOut:
        if sentFile is None:
            sys.stderr.write("Error: .gz output enabled, but output filehandle is null!\n")
            return
        else:
            if original:
                sentFile.write("[%s] ||| %s\n"%(nodeLabel, phrase_pair))
            else: #print Hiero grammar
                if isLex and srcNTCounter < 3:
                    nodeStr = "[S]" if nodeLabel == 0 else "[X]"
                    sentFile.write("[%s] ||| %s\n"%(nodeStr, phrase_pair))
                elif srcNTCounter > 2:
                    sys.stderr.write("Sentence %d:Rule '[X] ||| %s' has more than 2 NTs\n"%(count, phrase_pair))
    else: #printing to stdout
        if original:
            print "[%s] ||| %s"%(nodeLabel, phrase_pair)
        else:
            if isLex and srcNTCounter < 3:
                nodeStr = "S" if nodeLabel == 0 else "X"
                print "[%s] ||| %s"%(nodeStr, phrase_pair)
            elif srcNTCounter > 2:
                sys.stderr.write("Sentence %d:Rule '[X] ||| %s' has more than 2 NTs\n"%(count, phrase_pair))

def main():
    global NTLabel
    (opts, args) = getopt.getopt(sys.argv[1:], 'dz:')
    original = False
    zipOut = False
    dirLoc = ""
    for opt in opts: #flag/option handling
        if opt[0] == '-d': #prints out full output
            original = True
        elif opt[0] == '-z': #prints out in .gz
            zipOut = True
            dirLoc = opt[1]
            if not os.path.exists(dirLoc): #if directory does not exist, create it
                os.makedirs(dirLoc) 
    src_tgt_lines = open(args[0], 'r').read().splitlines()
    src_sent, tgt_sent = zip(*(line.split(" ||| ") for line in src_tgt_lines)) #assuming input is formatted as source_sent ||| target_sent #separate lists for src and tgt
    alignment_lines = open(args[1], 'r').read().splitlines()
    alignment_dicts = createAlignmentDicts(alignment_lines)
    count = 0 #line counter
    for line in sys.stdin:
        NTLabel = 0 #reset NTLabel before parsing an s-expression
        src = src_sent[count].split() #tokenized into separate list elements
        tgt = tgt_sent[count].split()
        alignDict = alignment_dicts[count]
        sentFile = None
        if (zipOut): #set up the output if need be
            sentLoc = dirLoc + "/grammar.%d.gz"%(count) #if we need to write out to zipped grammar, open file
            sentFile = gzip.open(sentLoc, 'wb')
        extractRules(count, original, line.strip(), src, tgt, alignDict, NTLabel, zipOut, sentFile)
        count += 1
        if (zipOut):
            sentFile.close()        

if __name__ == "__main__":
    main()
