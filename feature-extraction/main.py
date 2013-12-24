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
#Updated Usage: we do not write out feature vectors anymore, only feature names
#python main.py -a -r/usr0/home/avneesh/grammar-loc /usr0/home/avneesh/mingrammar-loc featNameOut rank paramOutFile
##########################

import sys, commands, string, gzip, os, os.path, re, getopt, cPickle
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import scipy.io as io
from tree import tree

exampleIDs = {}
exampleID = 0
inFeatIDs = {}
inFeatID = 0
outFeatIDs = {}
outFeatID = 0
def main():
    arityF = False
    realF = False
    hiero_loc = "" #for real-valued features
    lexF = False
    #any other features, put a boolean here
    (opts, args) = getopt.getopt(sys.argv[1:], 'ar:')
    for opt in opts:
        if opt[0] == '-a':
            arityF = True
        elif opt[0] == '-r':
            realF = True
            hiero_loc = opt[1]
        elif opt[0] == '-l':
            lexF = True
    minRule_grammars_loc = args[0]
    numSentences = len(os.listdir(minRule_grammars_loc))
    featname_out_loc = args[1]
    inFeatures = []
    outFeatures = []
    root_rules = [] #maintain row indices of rules that occur at root
    for line_num in range(0, numSentences): #loop through all sentence pairs
        minrule_fh = gzip.open(minRule_grammars_loc + "grammar.%d.gz"%(line_num))
        sync_tree = tree(0, None, None, minrule_fh) #use tree class to generate the tree using the minimal grammar for this sentence                        
        if realF: #if also extracting real-valued features, use the regular hiero grammar which has these features (from the suffix array)
            hiero_grammar_loc = hiero_loc + "/grammar.%d.gz"%(line_num) 
            root_rules.append(update_features(sync_tree, inFeatures, outFeatures, arityF, realF, hiero_grammar_loc))
        else:
            root_rules.append(update_features(sync_tree, inFeatures, outFeatures, arityF, realF))
    sys.stderr.write("Feature extraction complete\n")
    sys.stdout.flush()
    kappa = 5.0
    rank = int(args[3])
    in_fm = convertToSpMat(inFeatures, len(inFeatIDs), kappa) #convert to SpMat also does the feature scaling
    out_fm = convertToSpMat(outFeatures, len(outFeatIDs), kappa)            
    Y,Z = SVDandProjection(in_fm, out_fm, rank, True) #compute avg outer product, call Matlab to do SVD, extract Y and Z matrices
    sys.stderr.write("SVD and projections complete\n")    
    paramDict = computeCorrelations(Y,Z) #compute tensors, matrices, and vectors
    paramDict['Pi'] = estimatePiParams(root_rules, Y, rank) 
    sys.stderr.write("Parameter estimation complete\n")
    cPickle.dump(paramDict, open(args[2], "wb")) #write out in cPickle format for I/O algorithm
    featNameHandle = open(featname_out_loc, 'w')
    for feature in inFeatIDs: #write out feature names and IDs
        print >> featNameHandle, "INSIDE %s:%d"%(feature, inFeatIDs[feature])
    for feature in outFeatIDs:
        print >> featNameHandle, "OUTSIDE %s:%d"%(feature, outFeatIDs[feature])
    featNameHandle.close()     

def printTree(stree):
    print "Current rule: "
    print stree.rule
    for child in stree.children:
        print "going to child: "
        print child
        printTree(child)
        print "back from child: "
        print child

def SVDandProjection(inFeatMat, outFeatMat, rank, matlab):
    avgOP = (1.0 / inFeatMat.shape[0]) * (inFeatMat.transpose() * outFeatMat)
    U, S, V = matlabInterface(avgOP, rank) if matlab else la.svd(avgOP, full_matrices=False)    
    Y = inFeatMat * U #inFeatMat is a sparse matrix, so * is overloaded to represent matrix mult!!
    Z = outFeatMat.dot(V).dot(np.linalg.inv(S))
    return (Y, Z)

def matlabInterface(avgOP, rank):
    pwd = os.getcwd()
    out_loc = pwd + "/matlab_temp"
    io.savemat(out_loc, {'avgOP': avgOP})
    path = os.path.abspath(os.path.dirname(sys.argv[0]))
    os.chdir(path)
    os.system('matlab -nodesktop -nosplash -nojvm -r "matlab_svd ' + out_loc + " %s"%rank + '"')
    os.chdir(pwd)
    mat_return = io.loadmat(out_loc)
    #os.remove(out_loc)
    return mat_return['U'].newbyteorder('='), mat_return['S'].newbyteorder('='), mat_return['V'].newbyteorder('=')

def estimatePiParams(root_rules, Y, rank):
    outerProd = np.zeros(shape=(rank))
    for row_idx in root_rules:
        outerProd += Y[row_idx,:].transpose()
    outerProd = np.multiply(outerProd, 1.0/len(root_rules))
    return outerProd

def computeCorrelations(Y, Z):    
    paramDict = {}
    numExamples = Y.shape[0]
    rank = Y.shape[1]
    for rule in exampleIDs: #rule is an actual rule LHS ||| src ||| tgt
        arity = len(exampleIDs[rule][0].split())
        scale = 1.0 / numExamples
        outerProd = None
        if arity == 3: #compute tensor
            outerProd = np.zeros(shape=(rank, rank, rank))
            for ins_out_pair in exampleIDs[rule]:
                in_idx = [int(num) for num in re.findall("In:(\d+)", ins_out_pair)]
                out_idx = int(re.findall("Out:(\d+)", ins_out_pair)[0])
                outerProd += tensorProduct(Z[out_idx,:], Y[in_idx[0],:], Y[in_idx[1],:])            
        elif arity == 2: #compute matrix
            outerProd = np.zeros(shape=(rank, rank))
            for ins_out_pair in exampleIDs[rule]:
                in_idx = int(re.findall("In:(\d+)", ins_out_pair)[0])
                out_idx = int(re.findall("Out:(\d+)", ins_out_pair)[0])
                outerProd += np.outer(Z[out_idx,:], Y[in_idx,:])
        elif arity == 1:
            outerProd = np.zeros(shape=(rank))
            for ins_out_pair in exampleIDs[rule]:
                out_idx = int(re.findall("Out:(\d+)", ins_out_pair)[0])
                outerProd += Z[out_idx,:]
        else:
            sys.stderr.write("Rule has more than 2 inside trees or more than 1 outside tree!\n%s\n"%rule)
        outerProd = np.multiply(outerProd, scale) #scale by MLE count #is outerprod a matrix or array??
        paramDict[rule] = outerProd                       
    return paramDict
        
def tensorProduct(vec1, vec2, vec3):
    vs = [vec1, vec2, vec3]
    shape = map(len, vs)
    #specify the orientation of each vector before multiplication
    newshapes = np.diag(np.array(shape)-1)+1
    reshaped = [x.reshape(y) for x,y in zip(vs, newshapes)]
    #then, multiply directly
    return reduce(np.multiply, np.ix_(*vs))

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
    scaledFeat = rescaleFeatures(sparseFeat, kappa)
    return scaledFeat

def rescaleFeatures(sparseFeat, kappa):
    spFeatTransSquared = sparseFeat.transpose(copy=True)
    spFeatTransSquared.data **= 2
    scaleDenom = spFeatTransSquared.sum(axis=1) + kappa
    scaleVec = np.sqrt((sparseFeat.shape[0] - 1) * np.reciprocal(scaleDenom))
    scaleVecSp = sp.spdiags(scaleVec.flatten(), [0], len(scaleVec), len(scaleVec))
    return sparseFeat * scaleVecSp

def update_features(sent_tree, inFeat, outFeat, arityF, realF, hiero_loc=""):
    global exampleID, exampleIDs
    curID = 0
    if len(sent_tree.children) < 3: #filter rules with # of NTs > 2
        #check_maxLex(sent_tree)
        insFeatDict = {} #keys are feature IDs, vals are binary 1/0
        outFeatDict = {}
        extractRuleFeatures(sent_tree, insFeatDict, outFeatDict)
        if arityF:
            extractArityFeatures(sent_tree, insFeatDict, outFeatDict)
        if realF:
            extractRealValFeatures(hiero_loc, sent_tree, insFeatDict, outFeatDict)
        #any additional feature, define here
        inFeat.append(insFeatDict)
        outFeat.append(outFeatDict)
        curID = exampleID
        exampleID += 1
    exIDs = []
    for child in sent_tree.children: #at this stage, recurse on children
        childID = update_features(child, inFeat, outFeat, arityF, realF, hiero_loc)
        exIDs.append("In:%d"%(childID))
    if len(sent_tree.children) < 3:
        exIDs.append("Out:%d"%curID)
        key = sent_tree.rule #key is now the entire rule; to distinguish between S --> X1 X2 and X --> X1 X2
        #key = (sent_tree.src, sent_tree.tgt) #key is src-tgt phrase pair OLD!
        exampleList = exampleIDs[key] if key in exampleIDs else []
        exampleList.append(' '.join(exIDs))
        exampleIDs[key] = exampleList #link src/tgt phrase pair with the inside and outside row IDs
        print "%d %s ||| %s"%(curID, sent_tree.rule, ' '.join(exIDs))
    return curID

def check_maxLex(sent_tree):
    global max_Lex
    if len(sent_tree.children) < 1: #pre-terminal
        lenLex = len(sent_tree.src.split())
        if lenLex > max_Lex:
            max_Lex = lenLex
    else:
        plainLex = re.sub(r' \[.*?\] ', ' ', sent_tree.src) #remove NTs
        lenLex = len(plainLex.split())
        if lenLex > max_Lex:
            max_Lex = lenLex

def printFeatureVector(filehandle, featureDict, ruleID,  whichtree):
    featureVec = ["%d:%s"%(feature,str(featureDict[feature])) for feature in featureDict]
    print >> filehandle, "%d %s X %s"%(ruleID, whichtree, ' '.join(featureVec))

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
