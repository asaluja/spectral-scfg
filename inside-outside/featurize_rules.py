#!/usr/bin/python -tt

'''
File: featurize_rules.py 
Date: January 4, 2013 (based on decorate_minrules_with_features.py script)
Author: Avneesh Saluja (avneesh@cs.cmu.edu)
Description: this script takes as arguments the following:
arg 1: location of per-sentence grammars of minimal rules with associated spectral features
arg 2: location of files with full set of rules extracted from suffix array structure
(these rules are already featurized)
arg 3: location of output directory for featurized per-sentence grammars
arg 4: number of partitions
arg 5: this particular partition number
Usage: python featurize_rules.py minrules-dir fullrules-dir output-dir numPartitions partitionNum
'''

import sys, commands, string, gzip, os, re
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

def decorateSentenceGrammar(minrule_file, hiero_file, out_file, featureStr):
    hiero_fh = gzip.open(hiero_file, 'rb')
    hiero_rules = process_hiero_rules(hiero_fh)
    out_fh = gzip.open(out_file, 'w')
    if os.path.isfile(minrule_file):
        minrule_fh = gzip.open(minrule_file, 'rb')
        numRulesInHiero = 0
        numRulesTotal = 0
        for rule in minrule_fh:
            numRulesTotal += 1
            elements = rule.strip().split(' ||| ')
            key, noLex = removeSpanInfo(elements[:3])
            if key in hiero_rules:
                numRulesInHiero += 1
            ruleToPrint = rule.strip()
            if elements[1] == "<unk>": #src is OOV
                #ruleToPrint = "%s ||| %s ||| %s ||| %s OOV=1"%(elements[0], elements[2], elements[2], elements[3])
                ruleToPrint = "%s ||| %s ||| %s ||| %s"%(elements[0], elements[2], elements[2], elements[3])
            line_out = ruleToPrint + " DeletionRule=0.0 " + hiero_rules[key] if key in hiero_rules else ruleToPrint + " DeletionRule=1.0 " + featureStr
            if noLex: 
                assert key not in hiero_rules #if not lexical, key should never be in hiero rules
                line_out += " Glue=1.0"
                ntNumbers = [int(ntIdx) for ntIdx in re.findall(r'\[([^]]*)\]', elements[2])]
                if len(ntNumbers) == 2 and (ntNumbers[0] > ntNumbers[1]):                        
                    line_out += "Inverse=1"
            if elements[0] == '[S]':
                out_fh.write("%s ||| 0\n"%(' ||| '.join(elements[:3])))
            else:
                out_fh.write("%s\n"%line_out)
        minrule_fh.close()
        print "for minrule file %s, out of %d rules, %d are also in hiero"%(minrule_file, numRulesTotal, numRulesInHiero)
    else:
        print "could not find minrule grammar for sentence number %d, taking hiero grammar instead"%(sentNum)
        for key in hiero_rules:
            out_fh.write("%s ||| spectral=0.0 DeletionRule=0.0 %s\n"%(key, hiero_rules[key]))
    hiero_fh.close()
    out_fh.close()
        

def main():
    args = sys.argv[1:]
    minrule_grammars_loc = args[0]
    hiero_grammars_loc = args[1]
    outDir_loc = args[2]    
    numProcesses = int(args[3])
    minrule_grammars = os.listdir(minrule_grammars_loc)
    hiero_grammars = os.listdir(hiero_grammars_loc)
    features = extractFeatureList(hiero_grammars_loc + hiero_grammars[0])
    featureStr = ' '.join(["%s=0.0"%(feature) for feature in features])
    pool = mp.Pool(numProcesses)
    for minrule_file in minrule_grammars:
        pool.apply_async(decorateSentenceGrammar, (minrule_grammars_loc + minrule_file, hiero_grammars_loc + minrule_file, outDir_loc + minrule_file, featureStr))
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
