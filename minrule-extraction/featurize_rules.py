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
Update (Jan 4, 2014): included per-sentence grammar writing ability
Update (Jan 6, 2014): included filtering capability based on MLE for non per-sentence grammar writing
Update (Jan 7, 2014): made a multiprocess version of this code
Update (May 16, 2014): combined this script with inside-outside/featurize_rules.py, and renamed it
featurized_rulespy
Usage: python featurize_rules.py minrules-dir/spectral-marginal-dir hiero-dir/ output-dir/ numPartitions 
There are also 3 optional flags:
-f N: filter rules by P(e|f).  The argument N is how many rules to keep for a given source RHS. -f and -m
cannot be on together. 
-m: meaning the input grammar is the output of intersect_scfg.py (contains marginals), and not the output
of tree_to_rule.py (tool to convert ZGC minimal grammar extractor output to input of feature_extraction.py
or a format that cdec can read in).  
-s: write out a per sentence grammar.  This is done by default if -m is one.  -f and -s cannot be on 
together, because we need to aggregate information. 
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

def removeSpanInfo(elements):
    noLex = False
    expr = re.compile(r'\[([^]]*)\]')
    LHS = elements[0] if elements[0] == '[S]' else elements[0].split('_')[0] + ']'
    srcRHSArr = []    
    NTCount = 1
    for item in elements[1].split():
        if expr.match(item):
            srcRHSArr.append("[X,%d]"%NTCount)
            NTCount += 1
        else:
            srcRHSArr.append(item)
    if NTCount - 1 == len(elements[1].split()): #then we are dealing with a rule with no lex items
        noLex = True
    tgtRHSArr = []
    for item in elements[2].split():
        if expr.match(item):
            tgtRHSArr.append("[X," + item[1:])
        else:
            tgtRHSArr.append(item)
    srcRHS = ' '.join(srcRHSArr)    
    tgtRHS = ' '.join(tgtRHSArr) 
    if srcRHS == "<unk>": #OOV
        srcRHS = tgtRHS #tgtRHS has the source word, we passed this through in hg_io.py
    return (' ||| '.join([LHS, srcRHS, tgtRHS]), noLex)

def filterRules(countDict, limit):
    srcTgtDict = {}
    for key in countDict.keys(): #first, need to reprocess countDict
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

def computeFeatures(ruleToPrint):
    x = 2

def decorateSentenceGrammar(minRule_file, hiero_file, out_file, featureStr, optDict):
        hiero_fh = gzip.open(hiero_file, 'rb')
        hiero_rules = process_hiero_rules(hiero_fh)
        out_fh = None
        perSent = "perSentence" in optDict
        marginal = "marginal" in optDict
        if perSent or marginal:
            out_fh = gzip.open(out_file, 'w')
        numRulesInHiero = 0
        numRulesTotal = 0
        if os.path.isfile(minRule_file):
            minrule_fh = gzip.open(minRule_file, 'rb')
            for rule in minrule_fh:
                numRulesTotal += 1
                elements = rule.strip().split(' ||| ')
                noLex = False
                key = ' ||| '.join(elements[:3])
                if marginal:
                    key, noLex = removeSpanInfo(elements[:3])
                ruleToPrint = rule.strip() #assumption is that rule will be of the form a ||| b ||| c |||
                if elements[1] == "<unk>":
                    ruleToPrint = "%s ||| %s ||| %s ||| %s PassThrough=1"%(elements[0], elements[2], elements[2], elements[3])           
                if key in hiero_rules:
                    numRulesInHiero += 1
                    ruleToPrint += hiero_rules[key]
                else:
                    ruleToPrint += " minRule=1.0 "
                    if not noLex: #if featurizing without marginals, then we only read in noLex items so its true; otherwise, determined above
                        ruleToPrint += computeFeatures(ruleToPrint)
                if marginal and noLex:
                    assert key not in hiero_rules #if not lexical, key should never be in hiero rules
                    line_out += " Glue=1.0"
                    ntNumbers = [int(ntIdx) for ntIdx in re.findall(r'\[([^]]*)\]', elements[2])]
                    if len(ntNumbers) == 2 and (ntNumbers[0] > ntNumbers[1]):                        
                        line_out += " Inverse=1.0"
                if marginal:
                    if elements[0] == '[S]': #reached the end, so we should close? 
                        out_fh.write("%s ||| 0\n"%(' ||| '.join(elements[:3])))
                        out_fh.close()
                    else:
                        out_fh.write("%s\n"%line_out)
                else:
                    if "perSentence" in optDict:
                        out_fh.write("%s\n"%line_out)
                    else: #maintain a count of each src_tgt rule
                        seen_rules.append(line_out)
                        srcTgtCount = countDict[key] if key in countDict else 0
                        countDict[key] = srcTgtCount + 1
            minrule_fh.close()
            print "for grammar %s, out of %d rules, %d are also in hiero"%(minRule_file, numRulesTotal, numRulesInHiero)
        else:
            print "could not find minrule grammar for sentence number %d, taking hiero grammar instead"%(sentNum)
            for key in hiero_rules:
                out_fh.write("%s ||| spectral=0.0 DeletionRule=0.0 %s\n"%(key, hiero_rules[key]))
        hiero_fh.close()
        if perSent and not marginal: #add the NT only rules, but only if we're not reading in marginals (since we already have the NT rules), then close
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
    (opts, args) = getopt.getopt(sys.argv[1:], 'f:ms')
    for opt in opts:
        if opt[0] == '-f':
            optDict["filterRules"] = int(opt[1])
        elif opt[0] == '-m': #output of intersect_scfg.py
            optDict["marginal"] = 1
        elif opt[0] == '-s': #write out per sentence
            optDict["perSentence"] = 1
    if "filterRules" in optDict and ("perSentence" in optDict or "marginal" in optDict):
        sys.stderr.write("Error: -f and -s or -f and -m cannot be on at the same time\n")
        sys.exit()
    minRule_grammars_loc = args[0]
    hiero_grammars_loc = args[1]
    outFile_loc = args[2]
    numProcesses = int(args[3])
    minRule_grammars = os.listdir(minRule_grammars_loc)
    hiero_grammars = os.listdir(hiero_grammars_loc)        
    features = extractFeatureList(hiero_grammars_loc + hiero_grammars[0])
    featureStr = ' '.join(["%s=0.0"%(feature) for feature in features])
    seen_rules = None
    countDict = None
    pool = None
    if "marginal" not in optDict:
        seen_rules = mp.Manager().list()
        countDict = mp.Manager().dict()
        pool = mp.Pool(processes=numProcesses, initializer=init, initargs=(seen_rules, countDict))
    else:
        pool = mp.Pool(numProcesses)
    for minRule_file in minRule_grammars:
        pool.apply_async(decorateSentenceGrammar, (minRule_grammars_loc + minRule_file, hiero_grammars_loc + minRule_file, outFile_loc + minRule_file, featureStr, optDict))    
    pool.close()
    pool.join()

    if "marginal" not in optDict:
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
            output_fh.write("[X] ||| [X] [X] ||| [1] [2] ||| DeletionRule=1.0 Glue=1 %s\n"%featureStr)
            output_fh.write("[X] ||| [X] [X] ||| [2] [1] ||| DeletionRule=1.0 Glue=1 Inverse=1 %s\n"%featureStr)        
            output_fh.write("[S] ||| [X] ||| [1] ||| 0\n") #no features defined on the top-level rule, just for parsing completion purposes
            output_fh.close()

if __name__ == "__main__":
    main()
