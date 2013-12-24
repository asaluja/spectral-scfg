#!/usr/bin/python -tt

##########################
#File: main.py
#Description: main file for feature extraction
#arg1: directory location of minimal grammar files
#the files should be generated from the "tree_to_rule.py" script with the -d and -z flags
#arg2: file location of feature vectors for the rules
#stdout: rules (in the hiero format) decorated with pointers to the inside and outside feature vectors
#Date: October 16, 2013
#Update: November 4, 2013: added options to script, and new argument
#arg3: location of featureIDs and feature names to be written out
#In addition, we have added options corresponding to real-valued and arity features
#Usage: python main.py -a -r/usr0/home/avneesh/grammar-loc /usr0/home/avneesh/mingrammar-loc featVecOut featNameOut
##########################

import sys, commands, string, gzip, os, os.path, re, getopt
import numpy as np
import scipy.sparse as sp
import scipy.io as io
from tree import tree
from mlabwrap import mlab

ruleIDs = {}
ruleID = 0
inFeatIDs = {}
inFeatID = 0
outFeatIDs = {}
outFeatID = 0
def main():
    arityF = False
    realF = False
    #any other features, put a boolean here
    hiero_loc = "" #for real-valued features
    (opts, args) = getopt.getopt(sys.argv[1:], 'ar:')
    for opt in opts:
        if opt[0] == '-a':
            arityF = True
        elif opt[0] == '-r':
            realF = True
            hiero_loc = opt[1]
    minRule_grammars_loc = args[0]
    numSentences = len(os.listdir(minRule_grammars_loc))
    featname_out_loc = args[1]
    preT_OP_out = args[2]
    inT_OP_out = args[3]
    preT_in = []
    preT_out = []
    inT_in = []
    inT_out = []
    for line_num in range(0, numSentences): #loop through all sentence pairs
        minrule_fh = gzip.open(minRule_grammars_loc + "grammar.%d.gz"%(line_num))
        rule_dict = read_grammar_to_dict(minrule_fh) #grammar_dict is used to initialize sync_tree below
        sync_tree = tree('[A]', None, rule_dict)
        populate_ruleIDs(sync_tree)    
        if realF:
            hiero_grammar_loc = hiero_loc + "/grammar.%d.gz"%(line_num)
            update_features(sync_tree, preT_in, preT_out, inT_in, inT_out, arityF, realF, hiero_grammar_loc)
        else:
            update_features(sync_tree, preT_in, preT_out, inT_in, inT_out, arityF, realF)
    kappa = 5.0
    rank = 16
    preT_in_fm = convertToSpMat(preT_in, len(inFeatIDs), kappa)
    preT_out_fm = convertToSpMat(preT_out, len(outFeatIDs), kappa)
    inT_in_fm = convertToSpMat(inT_in, len(inFeatIDs), kappa)
    inT_out_fm = convertToSpMat(inT_out, len(outFeatIDs), kappa)    
    #printAvgOuterProduct(preT_in_fm, preT_out_fm, preT_OP_out)
    #printAvgOuterProduct(inT_in_fm, inT_out_fm, inT_OP_out)
    preT_result = computeSVD(preT_in_fm, preT_out_fm, rank)
    inT_result = computeSVD(inT_in_fm, inT_out_fm, rank)
    preT_features = projectFeatures(preT_in_fm, preT_out_fm, preT_result[0], preT_result[2])
    inT_features = projectFeatures(inT_in_fm, inT_out_fm, inT_result[0], inT_result[2])
    
    featNameHandle = open(featname_out_loc, 'w')
    for feature in inFeatIDs:
        print >> featNameHandle, "INSIDE %s:%d"%(feature, inFeatIDs[feature])
    for feature in outFeatIDs:
        print >> featNameHandle, "OUTSIDE %s:%d"%(feature, outFeatIDs[feature])
    featNameHandle.close()

def computeSVD(inFeatMat, outFeatMat, rank):
    avgOP = (1.0 / inFeatMat.shape[0]) * (inFeatMat.transpose() * outFeatMat)
    U, S, V = mlab.svds(avgOP, rank, nout=3)
    return (U, S, V)

def projectFeatures(in_fm, out_fm, U, V):
    Y = U.T * in_fm.T
    Z = V.T * out_fm.T
    return (Y.T, Z.T)

def computeCorrelations():
    x = 2

def printAvgOuterProduct(inFeatMat, outFeatMat, out_loc):
    avgOP = (1.0 / inFeatMat.shape[0]) * (inFeatMat.transpose() * outFeatMat)
    io.savemat(out_loc, {'avgOP': avgOP})

def convertToSpMat(rawFeatMat, numFeatures, kappa):
    rows = []
    cols = []
    vals = []
    for i, row in enumerate(rawFeatMat):
        rows.extend([i] * len(row))
        col, val = zip(*row.items())
        cols.extend(col)
        vals.extend(val)    
    sparseFeat = sp.csc_matrix((vals, (rows, cols)), shape=(len(rawFeatMat), numFeatures))
    return rescaleFeatures(sparseFeat, kappa)

def rescaleFeatures(sparseFeat, kappa):
    spFeatTransSquared = sparseFeat.transpose(copy=True)
    spFeatTransSquared.data **= 2
    scaleDenom = spFeatTransSquared.sum(axis=1) + kappa
    scaleVec = np.sqrt((sparseFeat.shape[0] - 1) * np.reciprocal(scaleDenom))
    scaleVecSp = sp.spdiags(scaleVec.flatten(), [0], len(scaleVec), len(scaleVec))
    return sparseFeat * scaleVecSp

def printFeatureVector(filehandle, featureDict, ruleID,  whichtree):
    featureVec = ["%d:%s"%(feature,str(featureDict[feature])) for feature in featureDict]
    print >> filehandle, "%d %s X %s"%(ruleID, whichtree, ' '.join(featureVec))

def read_grammar_to_dict(filehandle):
    rule_dict = {}
    for rule in filehandle:
        k = rule.strip().split(' ||| ')[0]
        raw_symbol = re.findall(r'\[([^]]?)\]', k)[0]
        if raw_symbol == "": #sometimes, the symbol is ']' which won't get picked up
            k = "[>]"
        rule_dict[k] = rule.strip()
    return rule_dict

def populate_ruleIDs(sent_tree):
    global ruleID, ruleIDs
    if sent_tree.rule not in ruleIDs:
        ruleIDs[sent_tree.rule] = ruleID
        ruleID += 1
    for child in sent_tree.children:
        populate_ruleIDs(child)

def update_features(sent_tree, preT_in, preT_out, inT_in, inT_out, arityF, realF, hiero_loc=""):
    #over here, we should probably filter rules out with numChildren > 2
    if len(sent_tree.children) < 3:
        insFeatDict = {}
        outFeatDict = {}
        extractRuleFeatures(sent_tree, insFeatDict, outFeatDict)
        if arityF:
            extractArityFeatures(sent_tree, insFeatDict, outFeatDict)
        if realF:
            extractRealValFeatures(hiero_loc, sent_tree, insFeatDict, outFeatDict)
        #any additional feature, define here
        if len(sent_tree.children) == 0: #pre-terminal
            preT_in.append(insFeatDict)
            preT_out.append(outFeatDict)
        else:
            inT_in.append(insFeatDict)
            inT_out.append(outFeatDict)
        printTrainingExample(sent_tree)
    for child in sent_tree.children: #at this stage, recurse on children
        update_features(child, preT_in, preT_out, inT_in, inT_out, arityF, realF, hiero_loc)

def printTrainingExample(sent_tree):
    insideVecIDs = []
    for child in sent_tree.children: #print out ruleID and corresponding inside and outside tree IDs
        insideVecIDs.append(ruleIDs[child.rule])        
    insideVecString = ' '.join(["In:%d"%(insideID) for insideID in insideVecIDs])
    outsideVecString = "" if sent_tree.parent is None else "Out:%d"%(ruleIDs[sent_tree.parent.rule])
    print "%d %s %s %s"%(ruleIDs[sent_tree.rule], sent_tree.rule, insideVecString, outsideVecString)
    
'''
Note: as in real-valued features, this is only defined on a particular
rule and not a tree.  So, this is only defined on the inside feature vector.
'''
def extractArityFeatures(sent_tree, inFeatDict, outFeatDict):
    numChildren = len(sent_tree.children)
    inFeatDict[addCheckFeature("arity_%d"%numChildren, "in")] = 1
    if numChildren > 1: #if so, check if the ordering is inverted
        if checkRuleInversion(sent_tree.rule):
            inFeatDict[addCheckFeature("reorder", "in")] = 1
    
def checkRuleInversion(rule):
    NT_numbers = [int(nt.split(',')[1]) for nt in re.findall(r'\[([^]]*)\]', rule.split(' ||| ')[2])]
    return False if sorted(NT_numbers) == NT_numbers else True

'''
Note: the real-valued features are associated with a particular rule, and
not with a tree.  Therefore, we choose to only define them for the inside
feature vector.  The outside feature dictionary is passed by convention, but ignored. 
'''
def extractRealValFeatures(hiero_loc, sent_tree, inFeatDict, outFeatDict):
    hiero_fh = gzip.open(hiero_loc, 'rb')
    hiero_rules = process_hiero_rules(hiero_fh)
    if sent_tree.rule in hiero_rules: #hiero rules only contain lexical item rules, so this means rule is lex
        rv_features = hiero_rules[sent_tree.rule] #I have a feeling this is the slow line
        for feature in rv_features.split():
            featElements = feature.split('=')
            featName = featElements[0]
            featVal = float(featElements[1])
            if featVal > 0:
                inFeatDict[addCheckFeature(featName, "in")] = featVal
        
def process_hiero_rules(filehandle):
    hiero_rules = {}
    for rule in filehandle:
        elements = rule.strip().split(' ||| ')
        key = ' ||| '.join(elements[:3])
        hiero_rules[key] = elements[3]
    return hiero_rules

'''
Note: the inside rule feature below extract the rule ID of the current rule,
as well as the rule of any children if they exist (distinguishing between
left and right).  The outside rule feature looks at the rule ID of the parent only.
'''
def extractRuleFeatures(sent_tree, inFeatDict, outFeatDict):
    inFeatDict[addCheckFeature("RuleSelf_%s"%sent_tree.rule, "in")] = 1 #self-rule ID feature
    child_num = 0
    for child in sent_tree.children: #then add children's IDs
        rule_format = "RuleLeft_%s"%(child.rule) if child_num == 0 else "RuleRight_%s"%(child.rule)
        child_num += 1
        inFeatDict[addCheckFeature(rule_format, "in")] = 1
    rule_format = "RuleParent_Root" if sent_tree.parent is None else "RuleParent_%s"%(sent_tree.parent.rule)
    outFeatDict[addCheckFeature(rule_format, "out")] = 1

def addCheckFeature(featStr, inOrOut):
    global inFeatID, inFeatIDs, outFeatID, outFeatIDs
    if inOrOut == "in":
        if featStr not in inFeatIDs:
            inFeatIDs[featStr] = inFeatID
            inFeatID += 1
        return inFeatIDs[featStr]
    else:
        if featStr not in outFeatIDs:
            outFeatIDs[featStr] = outFeatID
            outFeatID += 1
        return outFeatIDs[featStr]
        

if __name__ == "__main__":
    main()
