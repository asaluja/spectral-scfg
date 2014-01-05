#!/usr/bin/python -tt

'''
File: decorate_minrules_with_features.py
Date: October 13, 2013
Author: Avneesh Saluja (avneesh@cs.cmu.edu)
Description: this script takes as arguments the following:
arg 1: location of files with rules generated through the post-processed minimal
rule output (of of H. Zhang, D. Gildea, and D. Chiang, NAACL 2008)
arg 2: location of files with full set of rules extracted from suffix array structure
(these rules are already featurized)
arg 3: location of output files 
Usage: python decorate_minrules_with_features.py minRules-dir fullRules-dir output-dir
'''

import sys, commands, string, gzip, os, re

def process_hiero_rules(filehandle):
    hiero_rules = {}
    for rule in filehandle:
        elements = rule.strip().split(' ||| ')
        key = ' ||| '.join(elements[:3])
        hiero_rules[key] = elements[3]
    return hiero_rules

def extractFeatureList(gzfile):
    file_contents = gzip.open(gzfile, 'rb').readlines()
    features_line = file_contents[0].strip().split(' ||| ')[3]    
    features_list = [featPair.split('=')[0] for featPair in features_line.split()]
    return features_list

def removeSpanInfo(elements):
    expr = re.compile(r'\[([^]]*)\]')
    LHS = elements[0].split('_')[0] + ']'
    srcRHSArr = []
    NTCount = 1
    for item in elements[1].split():
        if expr.match(item):
            srcRHSArr.append("[X,%d]"%NTCount)
            NTCount += 1
        else:
            srcRHSArr.append(item)
    tgtRHSArr = []
    for item in elements[2].split():
        if expr.match(item):
            tgtRHSArr.append("[X," + item[1:])
        else:
            tgtRHSArr.append(item)
    srcRHS = ' '.join(srcRHSArr)
    tgtRHS = ' '.join(tgtRHSArr) 
    return ' ||| '.join([LHS, srcRHS, tgtRHS])

def main():
    args = sys.argv[1:]
    minRule_grammars_loc = args[0]
    hiero_grammars_loc = args[1]
    outDir_loc = args[2]
    numPartitions = int(args[3])
    partition = int(args[4])
    numSentences = len(os.listdir(minRule_grammars_loc))
    sentencesPerChunk = numSentences/numPartitions
    start = partition*sentencesPerChunk
    end = (partition+1)*sentencesPerChunk if partition < numPartitions - 1 else numSentences
    hiero_grammars = os.listdir(hiero_grammars_loc)        
    features = extractFeatureList(hiero_grammars_loc + hiero_grammars[0])
    featureStr = ' '.join(["%s=0.0"%(feature) for feature in features])
    for sentNum in range(start,end):
        hiero_fh = gzip.open(hiero_grammars_loc + "grammar.%d.gz"%sentNum, 'rb')
        hiero_rules = process_hiero_rules(hiero_fh)
        out_fh = gzip.open("%s/grammar.%d.gz"%(outDir_loc, sentNum), 'w')
        fileName = minRule_grammars_loc + "grammar.%d.gz"%sentNum
        if os.path.isfile(fileName):
            minrule_fh = gzip.open(fileName, 'rb')
            numRulesInHiero = 0
            numRulesTotal = 0
            for rule in minrule_fh:
                numRulesTotal += 1
                elements = rule.strip().split(' ||| ')
                key = removeSpanInfo(elements[:3])
                if key in hiero_rules:
                    numRulesInHiero += 1
                line_out = rule.strip() + " DeletionRule=0.0 " + hiero_rules[key] if key in hiero_rules else rule.strip() + " DeletionRule=1.0 " + featureStr
                out_fh.write("%s\n"%line_out)
            minrule_fh.close()
            print "for sentence number %d, out of %d rules, %d are also in hiero"%(sentNum, numRulesTotal, numRulesInHiero)
        else:
            for key in hiero_rules:
                out_fh.write("%s ||| spectral=0.0 DeletionRule=0.0 %s\n"%(key, hiero_rules[key]))
        out_fh.close()
        hiero_fh.close()

if __name__ == "__main__":
    main()
