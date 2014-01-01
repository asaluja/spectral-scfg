#!/usr/bin/python -tt

'''
File: feature_extraction.py
Author: Avneesh Saluja (avneesh@cs.cmu.edu)
Date: October 16, 2013
Description: main file for feature extraction
arg1: directory location of minimal grammar files
the files should be generated from the "tree_to_rule.py" script with the -d and -z flags
arg2: file location of feature vectors for the rules
stdout: rules (in the hiero format) decorated with pointers to the inside and outside feature vectors
Update: November 4, 2013: added options to script, and changed feature vector outputs in arg2 to featureIDs
arg2: location of featureIDs and feature names to be written out
In addition, we have added options corresponding to real-valued and arity features
By default, rule indicator features are always used. 
Usage: python feature_extraction.py -a -r/usr0/home/avneesh/feat.grammar-loc /usr0/home/avneesh/min.grammar-loc featNameOut rank paramOutFile
Update: December 31, 2013: added lexical and lexical class features.  Note that the word classes should 
be given in two separate files, one for each language, and that these paths should be separated by a ':'
Usage: python feature_extraction.py -a -c/path/to/wordClasses -r/path/to/suffixArrayGrammar /path/to/minimalGrammar /path/to/featureNamesOutput rank /path/to/parametersOutput
'''

import sys, commands, string, gzip, os, os.path, re, getopt, cPickle
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import scipy.io as io
from tree import tree

exampleIDs = {} #key: rule; values: list of strings, each string being a 1/2/3-tuple: InIdxLeft (if there), InIdxRight (if there), OutIdx
exampleID = 0
ruleToIdxMap = {} #key: rule; values: dict with two keys only: one key is the list of row indices corresponding to current rule that is key; the other is a list of row indices corresponding to parent
inFeatIDs = {} #key: feature (string), value: feature ID
inFeatID = 0
outFeatIDs = {} #same as inFeatIDs
outFeatID = 0
def main():
    featBinDict = {}
    (opts, args) = getopt.getopt(sys.argv[1:], 'ac:lr:')
    for opt in opts:
        if opt[0] == '-a':
            featBinDict["arityF"] = 1
        elif opt[0] == '-r':
            featBinDict["realF"] = opt[1]
        elif opt[0] == '-l':
            featBinDict["lexF"] = 1
        elif opt[0] == '-c':
            filenames = opt[1].split(':')
            featBinDict["classF"] = readInWordClasses(filenames)
    minRule_grammars_loc = args[0]
    numSentences = len(os.listdir(minRule_grammars_loc))
    featname_out_loc = args[1]
    inFeatures = []
    outFeatures = []
    root_rules = [] #maintain row indices of rules that occur at root
    for line_num in range(0, numSentences): #loop through all sentence pairs
        minrule_fh = gzip.open(minRule_grammars_loc + "grammar.%d.gz"%(line_num))
        sync_tree = tree(0, None, None, minrule_fh) #use tree class to generate the tree using the minimal grammar for this sentence                        
        root_rules.append(update_features(sync_tree, inFeatures, outFeatures, featBinDict))
    sys.stderr.write("Feature extraction complete\n")
    kappa = 5.0
    rank = int(args[2])    
    computeOOVProbMass(inFeatures)
    in_fm = convertToSpMat(inFeatures, len(inFeatIDs), kappa) #convert to SpMat also does the feature scaling
    out_fm = convertToSpMat(outFeatures, len(outFeatIDs), kappa)
    Y,Z = SVDandProjection(in_fm, out_fm, rank, True) #compute avg outer product, call Matlab to do SVD, extract Y and Z matrices
    sys.stderr.write("SVD and projections complete\n")    
    paramDict = computeCorrelations(Y,Z) #compute tensors, matrices, and vectors
    paramDict['Pi'] = estimatePiParams(root_rules, Y, rank) 
    sys.stderr.write("Parameter estimation complete\n")
    cPickle.dump(paramDict, open(args[3], "wb")) #write out in cPickle format for I/O algorithm
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

def readInWordClasses(filenames):
    srcClassDict = {}
    tgtClassDict = {}
    srcFile = filenames[0]
    tgtFile = filenames[1]
    for line in open(srcFile, 'rb'):
        elements = line.strip().split("\t")
        srcClassDict[elements[1]] = int(elements[0], 2)
    for line in open(tgtFile, 'rb'):
        elements = line.strip().split("\t")
        tgtClassDict[elements[1]] = int(elements[0], 2)
    return (srcClassDict, tgtClassDict)

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
        outerProd = np.multiply(outerProd, scale) #scale by MLE count
        elements = rule.split(' ||| ')
        src_key = ' ||| '.join(elements[:-1])
        tgt_key = elements[-1]
        srcDict = paramDict[src_key] if src_key in paramDict else {}
        srcDict[tgt_key] = outerProd
        paramDict[src_key] = srcDict
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

def computeOOVProbMass(inFeatures):
    global exampleIDs
    preTermSingletons = [k for k,v in exampleIDs.items() if len(v) == 1 and len(v[0].split()) == 1] #val is a string that provides current, inside and outside tree IDs
    print "Number of pre-term singletons: %d"%(len(preTermSingletons))
    numValidParents = 0
    for preTerm in preTermSingletons: #preTerm is a "LHS ||| src ||| tgt" string
        rowIdx = ruleToIdxMap[preTerm]["current"].pop()        
        featID = addCheckFeature("RuleSelf_%s"%preTerm, "in")        
        inFeatures[rowIdx].pop(featID)
        inFeatures[rowIdx][addCheckFeature("RuleSelf_OOV", "in")] = 1        
        if "parent" in ruleToIdxMap[preTerm]: #check if parent is included --> sometimes, parent is tertiary or higher rule which is invalid
            numValidParents += 1
            parentIdx = ruleToIdxMap[preTerm]["parent"].pop()        
            featIDLeft = addCheckFeature("RuleLeft_%s"%preTerm, "in")
            featIDRight = addCheckFeature("RuleRight_%s"%preTerm, "in")
            if featIDLeft in inFeatures[parentIdx]:
                inFeatures[parentIdx].pop(featIDLeft) #remove from dict
                inFeatures[parentIdx][addCheckFeature("RuleLeft_OOV", "in")] = 1
            elif featIDRight in inFeatures[parentIdx]:
                inFeatures[parentIdx].pop(featIDRight)
                inFeatures[parentIdx][addCheckFeature("RuleRight_OOV", "in")] = 1
            else:
                sys.stderr.write("Error! Cannot find parent of singleton preterminal %s\n"%preTerm)
        outIdxStr = exampleIDs[preTerm][0]
        exampleIDs.setdefault("[X] ||| <unk> ||| <unk>", []).append(outIdxStr)
        exampleIDs.pop(preTerm) #remove preTerm rule from exampleIDs
    print "Number of pre-term singletons with valid (unary or binary) parents: %d"%numValidParents

def update_features(sent_tree, inFeat, outFeat, featBinDict):
    global exampleID, exampleIDs
    global ruleToIdxMap
    curID = 0
    if len(sent_tree.children) < 3: #filter rules with # of NTs > 2
        insFeatDict = {} #keys are feature IDs, vals are binary 1/0
        outFeatDict = {}
        extractRuleFeatures(sent_tree, insFeatDict, outFeatDict)        
        if "arityF" in featBinDict:
            extractArityFeatures(sent_tree, insFeatDict, outFeatDict)
        if "realF" in featBinDict:
            extractRealValFeatures(featBinDict["realF"], sent_tree, insFeatDict, outFeatDict)
        if "lexF" in featBinDict:
            extractLexicalFeatures(sent_tree, insFeatDict, outFeatDict, False)
        if "classF" in featBinDict:
            extractLexicalFeatures(sent_tree, insFeatDict, outFeatDict, True, featBinDict["classF"])
        #define additional feature handling here
        inFeat.append(insFeatDict)
        outFeat.append(outFeatDict)
        curID = exampleID
        exampleID += 1
    exIDs = []
    for child in sent_tree.children: #at this stage, recurse on children
        childID = update_features(child, inFeat, outFeat, featBinDict)
        if len(child.children) == 0 and len(sent_tree.children) < 3: #means child rule is a preTerm
            ruleToIdxMap.setdefault(child.rule, {}).setdefault("parent", []).append(curID)
        exIDs.append("In:%d"%(childID))
    if len(sent_tree.children) < 3:
        exIDs.append("Out:%d"%curID)
        if len(sent_tree.children) == 0: #means current rule is a preTerm
            ruleToIdxMap.setdefault(sent_tree.rule, {}).setdefault("current", []).append(curID)
        key = sent_tree.rule #key is now the entire rule; to distinguish between S --> X1 X2 and X --> X1 X2
        exampleList = exampleIDs[key] if key in exampleIDs else []
        exampleList.append(' '.join(exIDs))
        exampleIDs[key] = exampleList #link src/tgt phrase pair with the inside and outside row IDs
        print "%d %s ||| %s"%(curID, sent_tree.rule, ' '.join(exIDs))
    return curID

def extractLexicalFeatures(sent_tree, inFeatDict, outFeatDict, isClass, classDictTuple=None):
    addRuleLexFeatures(sent_tree, inFeatDict, "in", isClass, classDictTuple) #incorporate isClass info here
    leftChild = True
    for child in sent_tree.children: #add yield features of children
        if isClass:
            srcY = [classDictTuple[0][word] for word in child.srcYield().split()] #all words in yield should be in classDict, because classDict is from training data
            tgtY = [classDictTuple[1][word] for word in child.tgtYield().split()]
            addYieldFeatures(srcY, tgtY, inFeatDict, leftChild, "in")
        else:
            addYieldFeatures(child.srcYield().split(), child.tgtYield().split(), inFeatDict, leftChild, "in")
        leftChild = False
    if sent_tree.parent is not None:
        addRuleLexFeatures(sent_tree.parent, outFeatDict, "out", isClass, classDictTuple)
        leftChild = True
        for child in sent_tree.parent.children:
            if child is not sent_tree: #the other child
                if isClass:
                    srcY = [classDictTuple[0][word] for word in child.srcYield().split()]
                    tgtY = [classDictTuple[1][word] for word in child.tgtYield().split()]
                    addYieldFeatures(srcY, tgtY, outFeatDict, leftChild, "out")
                else:
                    addYieldFeatures(child.srcYield().split(), child.tgtYield().split(), outFeatDict, leftChild, "out")
            else:
                leftChild = False

'''
Given a source and target yield in the form of lists, this function adds the first
and last words of the yields as features in the provided feature dictionary. 
'''
def addYieldFeatures(srcYield, tgtYield, featDict, isLeft, inOrOut):
    if isLeft:
        featDict[addCheckFeature("srcYieldLeftFirst_%s"%srcYield[0], inOrOut)] = 1
        featDict[addCheckFeature("srcYieldLeftLast_%s"%srcYield[-1], inOrOut)] = 1
        featDict[addCheckFeature("tgtYieldLeftFirst_%s"%tgtYield[0], inOrOut)] = 1
        featDict[addCheckFeature("tgtYieldLeftLast_%s"%tgtYield[-1], inOrOut)] = 1
    else:
        featDict[addCheckFeature("srcYieldRightFirst_%s"%srcYield[0], inOrOut)] = 1
        featDict[addCheckFeature("srcYieldRightLast_%s"%srcYield[-1], inOrOut)] = 1
        featDict[addCheckFeature("tgtYieldRightFirst_%s"%tgtYield[0], inOrOut)] = 1
        featDict[addCheckFeature("tgtYieldRightLast_%s"%tgtYield[-1], inOrOut)] = 1

'''
Given a rule, this function extracts the lexical items from the rule
(and if we are looking at word classes converts them into classes) and
adds the items as features for the given rule. 
'''
def addRuleLexFeatures(sent_tree, featDict, inOrOut, isClass, classDictTuple):
    expr = re.compile(r'\[([^]]*)\]')
    lexItemsSrc = [item for item in sent_tree.src.split() if not expr.match(item)]
    lexItemsTgt = [item for item in sent_tree.tgt.split() if not expr.match(item)]
    if isClass:
        lexItemsSrc = [classDictTuple[0][word] for word in lexItemsSrc]
        lexItemsTgt = [classDictTuple[1][word] for word in lexItemsTgt]
    for word in lexItemsSrc: #add lexical items of rule to inside features dict
        featDict[addCheckFeature("lexSrc_%s"%word, inOrOut)] = 1
    for word in lexItemsTgt:
        featDict[addCheckFeature("lexTgt_%s"%word, inOrOut)] = 1    
    
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
