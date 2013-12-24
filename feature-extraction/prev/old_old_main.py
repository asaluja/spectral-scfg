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
from tree import tree

ruleIDs = {}
ruleID = 0
featureIDs = {}
featureID = 0

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
    #feat_out_loc = args[1]
    featname_out_loc = args[1]
    preT_in = []
    preT_out = []
    inT_in = []
    inT_out = []
    #insFeatures = []
    #outFeatures = []        
    for line_num in range(0, numSentences): #loop through all sentence pairs
        minrule_fh = gzip.open(minRule_grammars_loc + "grammar.%d.gz"%(line_num))
        rule_dict = read_grammar_to_dict(minrule_fh)
        sync_tree = tree('[A]', None, rule_dict)
        #populate_ruleIDs(sync_tree, insFeatures, outFeatures) #add all new rules in this sentence
        populate_ruleIDs(sync_tree)    
        if realF:
            hiero_grammar_loc = hiero_loc + "/grammar.%d.gz"%(line_num)
            update_features(sync_tree, preT_in, preT_out, inT_in, inT_out, arityF, realF, hiero_grammar_loc)
            #update_features(sync_tree, insFeatures, outFeatures, arityF, realF, hiero_grammar_loc)
        else:
            #update_features(sync_tree, insFeatures, outFeatures, arityF, realF)
            update_features(sync_tree, preT_in, preT_out, inT_in, inT_out, arityF, realF)
    '''
    outhandle = open(feat_out_loc, 'w')
    for rule in ruleIDs: #loop through all ruleIDs
        featureDict = insFeatures[ruleIDs[rule]]
        printFeatureVector(outhandle, featureDict, ruleIDs[rule], "INSIDE")
        featureDict = outFeatures[ruleIDs[rule]]
        printFeatureVector(outhandle, featureDict, ruleIDs[rule], "OUTSIDE")
    outhandle.close()
    '''
    print preT_in
    print preT_out
    print inT_in
    print inT_out
    featNameHandle = open(featname_out_loc, 'w')
    for feature in featureIDs:
        print >> featNameHandle, "%s:%d"%(feature, featureIDs[feature])
    featNameHandle.close()

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

'''
def populate_ruleIDs(sent_tree, insFeatures, outFeatures):
    global ruleID, ruleIDs
    if sent_tree.rule not in ruleIDs:
        ruleIDs[sent_tree.rule] = ruleID
        ruleID += 1
        insFeatDict = {}
        outFeatDict = {}
        insFeatures.append(insFeatDict)
        outFeatures.append(outFeatDict)
    for child in sent_tree.children: #recursively populate rule IDs
        populate_ruleIDs(child, insFeatures, outFeatures)
'''

def update_features(sent_tree, preT_in, preT_out, inT_in, inT_out, arityF, realF, hiero_loc=""):
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
def update_features(sent_tree, insFeatures, outFeatures, arityF, realF, hiero_loc=""):    
    extractRuleFeatures(insFeatures[ruleIDs[sent_tree.rule]], outFeatures[ruleIDs[sent_tree.rule]], sent_tree)
    if realF:
        extractRealValFeatures(hiero_loc, insFeatures[ruleIDs[sent_tree.rule]], outFeatures[ruleIDs[sent_tree.rule]], sent_tree)
    if arityF:
        extractArityFeatures(insFeatures[ruleIDs[sent_tree.rule]], outFeatures[ruleIDs[sent_tree.rule]], sent_tree)
    #any additional feature, define here
    insideVecIDs = []
    for child in sent_tree.children: #print out ruleID and corresponding inside and outside tree IDs
        insideVecIDs.append(ruleIDs[child.rule])        
    insideVecString = ' '.join(["In:%d"%(insideID) for insideID in insideVecIDs])
    outsideVecString = "" if sent_tree.parent is None else "Out:%d"%(ruleIDs[sent_tree.parent.rule])
    print "%d %s %s %s"%(ruleIDs[sent_tree.rule], sent_tree.rule, insideVecString, outsideVecString)
    for child in sent_tree.children: #at this stage, recurse on children
        update_features(child, insFeatures, outFeatures, arityF, realF, hiero_loc)
'''

'''
Note: as in real-valued features, this is only defined on a particular
rule and not a tree.  So, this is only defined on the inside feature vector.
'''
def extractArityFeatures(sent_tree, inFeatDict, outFeatDict):
    numChildren = len(sent_tree.children)
    inFeatDict[addCheckFeature("arity_%d"%numChildren)] = 1
    if numChildren > 1: #if so, check if the ordering is inverted
        if checkRuleInversion(sent_tree.rule):
            inFeatDict[addCheckFeature("reorder")] = 1
    
'''
def extractArityFeatures(insFeatDict, outFeatDict, sent_tree):
    global featureID, featureIDs
    numChildren = len(sent_tree.children)
    featFormat = "arity_%d"%numChildren
    if featFormat not in featureIDs:
        featureIDs[featFormat] = featureID
        featureID += 1
    insFeatDict[featureIDs[featFormat]] = 1
    featFormat = "reorder"
    if featFormat not in featureIDs:
        featureIDs[featFormat] = featureID
        featureID += 1
    if numChildren > 1: #if so, check if the ordering is inverted
        if checkRuleInversion(sent_tree.rule):
            insFeatDict[featureIDs[featFormat]] = 1
'''

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
        for feature in features.split():
            featName = feature.split('=')[0]
            inFeatDict[addCheckFeature(featName)] = float(feature.split('=')[1])
'''
def extractRealValFeatures(hiero_loc, insFeatDict, outFeatDict, sent_tree):
    global featureID, featureIDs
    hiero_fh = gzip.open(hiero_loc, 'rb')
    hiero_rules = process_hiero_rules(hiero_fh)
    if sent_tree.rule in hiero_rules: #hiero rules only contain lexical item rules, so this means rule is lex
        features = hiero_rules[sent_tree.rule]
        for feature in features.split():
            featName = feature.split('=')[0]
            if featName not in featureIDs:
                featureIDs[featName] = featureID
                featureID += 1
            if not featureIDs[featName] in insFeatDict: 
                insFeatDict[featureIDs[featName]] = float(feature.split('=')[1])
'''
        
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
    inFeatDict[addCheckFeature("RuleSelf_%s"%sent_tree.rule)] = 1 #self-rule ID feature
    child_num = 0
    for child in sent_tree.children: #then add children's IDs
        rule_format = "RuleLeft_%s"%(child.rule) if child_num == 0 else "RuleRight_%s"%(child.rule)
        child_num += 1
        inFeatDict[addCheckFeature(rule_format)] = 1
    rule_format = "RuleParent_Root" if sent_tree.parent is None else "RuleParent_%s"%(sent_tree.parent.rule)
    outFeatDict[addCheckFeature(rule_format)] = 1

'''    
def extractRuleFeatures(insFeatDict, outFeatDict, sent_tree):
    global featureID, featureIDs
    if "RuleSelf_%s"%(sent_tree.rule) not in featureIDs: #first, add self rule identity    
        featureIDs["RuleSelf_%s"%sent_tree.rule] = featureID
        featureID += 1
    addCheckFeature(featureIDs["RuleSelf_%s"%sent_tree.rule], insFeatDict)
    child_num = 0
    for child in sent_tree.children: #then, add childrens' identities
        rule_format = "RuleLeft_%s"%(child.rule) if child_num == 0 else "RuleRight_%s"%(child.rule)
        child_num += 1
        if rule_format not in featureIDs:
            featureIDs[rule_format] = featureID
            featureID += 1
        addCheckFeature(featureIDs[rule_format], insFeatDict)        
    #now for outside features
    rule_format = "RuleParent_Root" if sent_tree.parent is None else "RuleParent_%s"%(sent_tree.parent.rule)
    if rule_format not in featureIDs:
        featureIDs[rule_format] = featureID
        featureID += 1
    addCheckFeature(featureIDs[rule_format], outFeatDict) 
'''

def addCheckFeature(featStr):
    global featureID, featureIDs
    if featStr not in featureIDs:
        featureIDs[featStr] = featureID
        featureID += 1
    return featureIDs[featStr]

'''
def addCheckFeature(featID, featDict):
    if featID in featDict:
        featDict[featID] += 1
    else:
        featDict[featID] = 1
'''        

if __name__ == "__main__":
    main()
