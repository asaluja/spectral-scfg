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
Update (Jan 4, 2013): included per-sentence grammar writing ability
Update (Jan 6, 2013): included filtering capability based on MLE for non per-sentence grammar writing
Update (Jan 7, 2013): made a multiprocess version of this code
'''

import sys, commands, string, gzip, os, getopt, re
import multiprocessing as mp

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

def filterRules(countDict, limit):
    #first, need to reprocess countDict
    srcTgtDict = {}
    for key in countDict.keys():
        elements = key.split(' ||| ')
        srcKey = ' ||| '.join(elements[:2])
        tgtKey = elements[2]
        srcCountDict = srcTgtDict[srcKey] if srcKey in srcTgtDict else {}
        srcCountDict[tgtKey] = countDict[key]
        srcTgtDict[srcKey] = srcCountDict
    for srcKey in srcTgtDict:
        numTgtRules = len(srcTgtDict[srcKey])
        if numTgtRules > limit: #then we need to filter rules
            sorted_tgtRules = sorted(srcTgtDict[srcKey], key=srcTgtDict[srcKey].get, reverse=True)
            rules_to_filter = sorted_tgtRules[limit:]
            for rule in rules_to_filter:
                srcTgtDict[srcKey].pop(rule)
            if len(rules_to_filter) > 0:
                sys.stderr.write("Source RHS: %s; out of %d rules, filtered %d\n"%(srcKey, numTgtRules, len(rules_to_filter)))
    return srcTgtDict

def decorateSentenceGrammar(minRule_file, hiero_file, out_file, featureStr, optDict):
        minrule_fh = gzip.open(minRule_file, 'rb')
        hiero_fh = gzip.open(hiero_file, 'rb')
        out_fh = None
        if "perSentence" in optDict:
            out_fh = gzip.open(out_file, 'w')
        hiero_rules = process_hiero_rules(hiero_fh)
        numRulesInHiero = 0
        numRulesTotal = 0
        for rule in minrule_fh:
            numRulesTotal += 1
            elements = rule.strip().split(' ||| ')
            key = ' ||| '.join(elements[:3])
            if key in hiero_rules:
                numRulesInHiero += 1
            line_out = rule.strip() + " ||| DeletionRule=0.0 " + hiero_rules[key] if key in hiero_rules else rule.strip() + " ||| DeletionRule=1.0 " + featureStr
            if "perSentence" in optDict:
                out_fh.write("%s\n"%line_out)
            else: #maintain a count of each src_tgt rule
                seen_rules.append(line_out)
                srcTgtCount = countDict[key] if key in countDict else 0
                countDict[key] = srcTgtCount + 1
        minrule_fh.close()
        hiero_fh.close()
        print "for grammar %s, out of %d rules, %d are also in hiero"%(minRule_file, numRulesTotal, numRulesInHiero)
        if "perSentence" in optDict: #add the NT only rules, then close
            out_fh.write("[X] ||| [X,1] [X,2] ||| [X,1] [X,2] ||| DeletionRule=1.0 Glue=1 %s\n"%featureStr)
            out_fh.write("[X] ||| [X,1] [X,2] ||| [X,2] [X,1] ||| DeletionRule=1.0 Glue=1 Inverse=1 %s\n"%featureStr)            
            out_fh.write("[S] ||| [X,1] ||| [X,1] ||| 0\n") #no features defined on the top-level rule, just for parsing completion purposes
            out_fh.close()

def init(sr, cd):
    global seen_rules, countDict
    seen_rules = sr
    countDict = cd

def main():
    optDict = {}
    (opts, args) = getopt.getopt(sys.argv[1:], 'f:s')
    for opt in opts:
        if opt[0] == '-f':
            optDict["filterRules"] = int(opt[1])
        elif opt[0] == '-s': #write out per sentence
            optDict["perSentence"] = 1
    if "filterRules" in optDict and "perSentence" in optDict: #cannot handle this scenario yet
        sys.stderr.write("Cannot filter rules and write per sentence grammar at the same time\n")
        sys.exit()
    minRule_grammars_loc = args[0]
    hiero_grammars_loc = args[1]
    outFile_loc = args[2]
    numProcesses = int(args[3])
    minRule_grammars = os.listdir(minRule_grammars_loc)
    hiero_grammars = os.listdir(hiero_grammars_loc)        
    features = extractFeatureList(hiero_grammars_loc + hiero_grammars[0])
    featureStr = ' '.join(["%s=0.0"%(feature) for feature in features])
    seen_rules = mp.Manager().list()
    countDict = mp.Manager().dict()
    pool = mp.Pool(processes=numProcesses, initializer=init, initargs=(seen_rules, countDict))
    for minRule_file in minRule_grammars:
        pool.apply_async(decorateSentenceGrammar, (minRule_grammars_loc + minRule_file, hiero_grammars_loc + minRule_file, outFile_loc + minRule_file, featureStr, optDict))    
    pool.close()
    pool.join()
    print "number of rules seen: %d"%len(seen_rules)
    if "perSentence" not in optDict:
        output_fh = gzip.open(outFile_loc, 'wb')
        seen_rules_uniq = list(set(seen_rules))
        if "filterRules" in optDict:
            filteredDict = filterRules(countDict, optDict["filterRules"])            
            for rule in seen_rules_uniq:
                srcKey = ' ||| '.join(rule.split(' ||| ')[:2])
                tgtKey = rule.split(' ||| ')[2]
                if tgtKey in filteredDict[srcKey]: #i.e., we haven't pruned it away
                    output_fh.write("%s\n"%(rule))
        else:
            for rule in seen_rules_uniq:
                output_fh.write("%s\n"%(rule))
        output_fh.write("[X] ||| [X,1] [X,2] ||| [X,1] [X,2] ||| DeletionRule=1.0 Glue=1 %s\n"%featureStr)
        output_fh.write("[X] ||| [X,1] [X,2] ||| [X,2] [X,1] ||| DeletionRule=1.0 Glue=1 Inverse=1 %s\n"%featureStr)        
        output_fh.write("[S] ||| [X,1] ||| [X,1] ||| 0\n") #no features defined on the top-level rule, just for parsing completion purposes
        output_fh.close()

if __name__ == "__main__":
    main()
