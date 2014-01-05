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

import sys, commands, string, gzip, os, getopt, re

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

def main():
    optDict = {}
    (opts, args) = getopt.getopt(sys.argv[1:], 'ds')
    for opt in opts:
        if opt[0] == '-s': #write out per sentence
            optDict["perSentence"] = 1
    minRule_grammars_loc = args[0]
    hiero_grammars_loc = args[1]
    outFile_loc = args[2]
    seen_rules = []
    minRule_grammars = len(os.listdir(minRule_grammars_loc))
    hiero_grammars = os.listdir(hiero_grammars_loc)        
    features = extractFeatureList(hiero_grammars_loc + hiero_grammars[0])
    featureStr = ' '.join(["%s=0.0"%(feature) for feature in features])
    for minRule_file in minRule_grammars:
        minrule_fh = gzip.open(minRule_grammars_loc + minRule_file, 'rb')
        hiero_fh = gzip.open(hiero_grammars_loc + minRule_file, 'rb')
        out_fh = None
        if "perSentence" in optDict:
            out_fh = gzip.open("%s/%s"%(outFile_loc, minRule_file), 'w')
        hiero_rules = process_hiero_rules(hiero_fh)
        numRulesInHiero = 0
        numRulesTotal = 0
        for rule in minrule_fh:
            numRulesTotal += 1
            elements = rule.strip().split(' ||| ')
            key = ' ||| '.join(elements[:3])
            if key in hiero_rules:
                numRulesInHiero += 1
            line_out = rule.strip() + " DeletionRule=0.0 " + hiero_rules[key] if key in hiero_rules else rule.strip() + " DeletionRule=1.0 " + featureStr
            if "perSentence" in optDict:
                out_fh.write("%s\n"%line_out)
            else:
                seen_rules.append(line_out)
        minrule_fh.close()
        hiero_fh.close()
        print "for sentence number %d, out of %d rules, %d are also in hiero"%(sentNum, numRulesTotal, numRulesInHiero)
        if "perSentence" in optDict:
            out_fh.close()
    if "perSentence" not in optDict:
        seen_rules_uniq = list(set(seen_rules))
        output_fh = gzip.open(outFile_loc, 'wb')
        for rule in seen_rules_uniq:
            output_fh.write("%s\n"%(rule))
        output_fh.close()

if __name__ == "__main__":
    main()
