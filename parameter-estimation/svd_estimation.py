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
Update: January 3, 2014: added computing OOV probability mass as an option
Update: April 14, 2014: added maximum likelihood estimates for parameter estimation as an option
Update: April 21, 2014: added parameter smoothing as an option 
Update: May 2, 2014: added span length feature as an option
Update: May 18, 2014: writes out raw count dictionary for MLE (so that EM can use this information for OOV)
'''

import sys, commands, string, gzip, os, os.path, re, getopt, cPickle, time
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
    (opts, args) = getopt.getopt(sys.argv[1:], 'ac:f:lLmor:s:')
    for opt in opts:
        if opt[0] == '-a': #arity feature
            featBinDict["arityF"] = 1
        elif opt[0] == '-c': #lexical class features
            filenames = opt[1].split(':')
            featBinDict["classF"] = readInWordClasses(filenames)
        elif opt[0] == '-f': #filter rules
            featBinDict["filterRules"] = int(opt[1])
        elif opt[0] == '-l': #lexical features
            featBinDict["lexF"] = 1
        elif opt[0] == '-L': #span length feature
            featBinDict["lengthF"] = 1
        elif opt[0] == '-m': #ML estimate
            featBinDict["MLE"] = opt[1] #argument writes out raw count dictionary before normalization
        elif opt[0] == '-o': #OOV
            featBinDict["OOV"] = 1
        elif opt[0] == '-r': #real-valued features
            featBinDict["realF"] = opt[1]
        elif opt[0] == '-s': #smoothing option
            featBinDict["smooth"] = float(opt[1])
    minRule_grammars_loc = args[0]
    numSentences = len(os.listdir(minRule_grammars_loc))
    featname_out_loc = args[1]

    inFeatures = []
    outFeatures = []
    root_rules = [] #maintain row indices of rules that occur at root
    start = time.clock()
    for line_num in range(0, numSentences): #loop through all sentence pairs and extract features
        minrule_fh = gzip.open(minRule_grammars_loc + "grammar.%d.gz"%(line_num))
        sync_tree = tree(0, None, None, minrule_fh) #use tree class to generate the tree using the minimal grammar for this sentence                        
        root_rules.append(update_features(sync_tree, inFeatures, outFeatures, featBinDict))
    timeTaken = time.clock() - start
    sys.stderr.write("Feature extraction complete. Time taken: %.3f sec\n"%timeTaken)
    kappa = 5.0
    rank = int(args[2])
    singletons = []
    if "OOV" in featBinDict: #only treat unigram OOVs; we re-organize inside features to reflect this
        singletons = computeOOVProbMass(inFeatures) #rename function
    in_fm = convertToSpMat(inFeatures, len(inFeatIDs), kappa) #convert to SpMat also does the feature scaling
    out_fm = convertToSpMat(outFeatures, len(outFeatIDs), kappa)
    start = time.clock()
    Y,Z = SVDandProjection(in_fm, out_fm, rank) #compute avg outer product, call Matlab to do SVD, extract Y and Z matrices
    timeTaken = time.clock() - start
    sys.stderr.write("SVD and projections complete. Time taken: %.3f sec\n"%timeTaken)    
    start = time.clock()
    paramDict, countDict = computeCorrelationsAndSmooth(Y, Z, featBinDict["smooth"]) if "smooth" in featBinDict else computeCorrelations(Y, Z) #compute tensors, matrices, and vectors

    #below code for conditional prob write out
    '''
    normalizer = 0
    for srcKey in paramDict:
        normalizer = sum([countDict[srcKey][tgtKey] for tgtKey in paramDict[srcKey]])
        for tgtKey in paramDict[srcKey]:
            cond_prob = countDict[srcKey][tgtKey] / float(normalizer)
            paramDict[srcKey][tgtKey] = np.multiply(paramDict[srcKey][tgtKey], cond_prob)
    '''
    #cond prob write out code end
    #below code for tensor debugging
    '''
    G = np.random.random_sample((rank, rank))
    G = G + G.transpose()
    Ginv = np.linalg.inv(G)
    for src_RHS in paramDict:
        if src_RHS == "Pi":
            paramDict[src_RHS] = np.absolute(paramDict[src_RHS]).dot(Ginv)
        else:
            for target_RHS in paramDict[src_RHS]:
                parameter = np.absolute(paramDict[src_RHS][target_RHS])
                arity = len(parameter.shape) - 1
                if arity == 0:
                    paramDict[src_RHS][target_RHS] = G.dot(parameter)
                elif arity == 1:
                    paramDict[src_RHS][target_RHS] = G.dot(parameter).dot(Ginv)
                elif arity == 2:
                    result = np.tensordot(G, parameter, axes=[1,0])
                    result = np.tensordot(Ginv, result, axes=[1,1]).swapaxes(0,1)
                    result = np.tensordot(Ginv, result, axes=[1,2]).swapaxes(1,2)
                    paramDict[src_RHS][target_RHS] = result
    '''    
    #tensor debugging end

    paramDict['Pi'] = estimatePiParams(root_rules, Y, rank) 
    assignProbMassToOOVs("OOV" in featBinDict, paramDict, countDict, singletons)
    timeTaken = time.clock() - start
    sys.stderr.write("Parameter estimation complete. Time taken: %.3f sec\n"%timeTaken),
    if "filterRules" in featBinDict:
        start = time.clock()
        filterRulesByCount(countDict, paramDict, featBinDict["filterRules"])        
        timeTaken = time.clock() - start
        sys.stderr.write("Filtered rules to keep top %d per source RHS.  Time taken: %.3f sec\n"%(featBinDict["filterRules"], timeTaken))
    if "MLE" in featBinDict: #if we are just dealing with MLE estimates, then write them out directly
        start = time.clock()
        normalizeCountDict(countDict)
        cPickle.dump(countDict, open(args[3], "wb"))
        timeTaken = time.clock() - start
        sys.stderr.write("Normalized counts, and wrote out normalized counts.  Time taken: %.3f sec\n"%timeTaken)
    else:
        start = time.clock()
        cPickle.dump(paramDict, open(args[3], "wb")) #write out in cPickle format for I/O algorithm
        featNameHandle = open(featname_out_loc, 'w')
        for feature in inFeatIDs: #write out feature names and IDs
            print >> featNameHandle, "INSIDE %s:%d"%(feature, inFeatIDs[feature])
        for feature in outFeatIDs:
            print >> featNameHandle, "OUTSIDE %s:%d"%(feature, outFeatIDs[feature])
        featNameHandle.close()        
        timeTaken = time.clock() - start
        sys.stderr.write("Wrote out parameters and features.  Time taken: %.3f sec\n"%timeTaken)
    sys.stderr.write("Parameter estimation complete\n")

'''
for debugging purposes
'''
def printTree(stree):
    print "Current rule: "
    print stree.rule
    for child in stree.children:
        print "going to child: "
        print child
        printTree(child)
        print "back from child: "
        print child

def assignProbMassToOOVs(OOV, paramDict, countDict, singletons):
    if not OOV: #assign 0 to OOVs if we explicitly do not want a parameter for them
        for singleton in singletons:
            elements = ' ||| '.split(singleton)
            srcKey = ' ||| '.join(elements[:2])
            tgtKey = elements[2]
            srcDict = paramDict[srcKey] if srcKey in paramDict else {}
            srcDict[tgtKey] = np.zeros((rank))
            paramDict[srcKey] = srcDict
            countSrcDict = countDict[srcKey] if srcKey in countDict else {}
            countSrcDict[tgtKey] = 0
            countDict[srcKey] = countSrcDict
    else:
        OOV_param = paramDict["[X] ||| <unk>"]["<unk>"]
        avgOOV_param = np.divide(OOV_param, len(singletons)+1)
        paramDict["[X] ||| <unk>"]["<unk>"] = avgOOV_param
        countDict["[X] ||| <unk>"]["<unk>"] = 1
        for singleton in singletons:
            elements = singleton.split(' ||| ')
            srcKey = ' ||| '.join(elements[:2])
            tgtKey = elements[2]
            srcDict = paramDict[srcKey] if srcKey in paramDict else {}
            srcDict[tgtKey] = avgOOV_param
            paramDict[srcKey] = srcDict
            countSrcDict = countDict[srcKey] if srcKey in countDict else {}
            countSrcDict[tgtKey] = 1
            countDict[srcKey] = countSrcDict

'''
for word class feature - reads in word classes that were
produced from P. Liang's Brown Clustering code (two files,
one for each language) and returns a pair of dictionaries
'''
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

'''
given two feature matrices and a rank, this code first computes
the average outer product of the features, then interfaces with matlab
to compute the SVD, and then returns lower-dimensional projections of
the feature matrices. 
'''
def SVDandProjection(inFeatMat, outFeatMat, rank):
    avgOP = (1.0 / inFeatMat.shape[0]) * (inFeatMat.transpose() * outFeatMat)
    U, S, V = matlabInterface(avgOP, rank)
    print "Singular values are: "
    print S.diagonal()
    Y = inFeatMat * U #inFeatMat is a sparse matrix, so * is overloaded to represent matrix mult!!
    Z = outFeatMat.dot(V).dot(np.linalg.inv(S))
    return (Y, Z)

'''
specific function using os system calls to interface
with matlab via the command line. 
'''
def matlabInterface(avgOP, rank):
    pwd = os.getcwd()
    out_loc = pwd + "/matlab_temp"
    io.savemat(out_loc, {'avgOP': avgOP})
    path = os.path.abspath(os.path.dirname(sys.argv[0]))
    os.chdir(path)
    os.system('matlab -nodesktop -nosplash -nojvm -r "matlab_svd ' + out_loc + " %s"%rank + '"')
    os.chdir(pwd)
    mat_return = io.loadmat(out_loc)
    return mat_return['U'].newbyteorder('='), mat_return['S'].newbyteorder('='), mat_return['V'].newbyteorder('=')

'''
function used to prune rules using MLE.  if the number of unique target
side rules for a given source side exceeds a limit, we sort the rules by
count and prune away the low count rules beyond the limit. 
'''
def filterRulesByCount(countDict, paramDict, limit):
    for src_rule in countDict:
        if len(countDict[src_rule]) > limit: #then we need to prune
            sorted_tgtRules = sorted(countDict[src_rule], key=countDict[src_rule].get, reverse=True)
            rules_to_filter = sorted_tgtRules[limit:]            
            original_count = len(countDict[src_rule])
            for rule in rules_to_filter:
                paramDict[src_rule].pop(rule) #remove it from the parameter srcDict
                countDict[src_rule].pop(rule) #remove it from the count srcDict
            if len(rules_to_filter) > 0: #N.B.: need to update this to output the correct number
                sys.stderr.write("Source RHS: %s; out of %d rules, filtered %d\n"%(src_rule, original_count, len(rules_to_filter)))

'''
estimate the start of sentence parameters; we loop
through all rules that occur at the root (i.e., left hand side
is [S]) and accumulate their inside scores.  
'''
def estimatePiParams(root_rules, Y, rank):
    outerProd = np.zeros(shape=(rank))
    for row_idx in root_rules:
        outerProd += Y[row_idx,:].transpose()
    outerProd = np.multiply(outerProd, 1.0/len(root_rules))
    #outerProd = np.multiply(outerProd, 1.0/Y.shape[0])
    return outerProd

'''
This function is exactly the same as computeCorrelations, except we also smooth
the parameters according to the smoothing scheme described in Cohen et al., NAACL 2013
but with suitable modifications to also handle unary rules. 
We do not do any smoothing for the lexical rules. 
'''
def computeCorrelationsAndSmooth(Y, Z, C):
    paramDict = {}
    countDict = {}
    numExamples = Y.shape[0]
    rank = Y.shape[1]
    F = np.multiply(np.sum(Y, axis=0), 1.0/numExamples)
    H = np.multiply(np.sum(Z, axis=0), 1.0/numExamples)
    E_4_tensor = tensorProduct(H, F, F)
    E_4_matrix = np.outer(H, F)
    for rule in exampleIDs: #rule is in the form 'LHS ||| src ||| tgt'
        arity = len(exampleIDs[rule][0].split())
        ruleCount = len(exampleIDs[rule]) #|Q^{a \rightarrow b c}| in the paper
        coeff = np.sqrt(ruleCount) / (C + np.sqrt(ruleCount))
        scale = float(ruleCount) / numExamples #MLE scaling for the end
        param, outerProd = None, None
        secMom_firstFree, secMom_secFree, secMom_thirdFree = None, None, None
        firstMom_firstFix, firstMom_secFix, firstMom_thirdFix = None, None, None
        if arity == 3: 
            outerProd = np.zeros(shape=(rank, rank, rank))
            secMom_firstFree, secMom_secFree, secMom_thirdFree = np.zeros(shape=(rank, rank)), np.zeros(shape=(rank, rank)), np.zeros(shape=(rank, rank))
            firstMom_firstFix, firstMom_secFix, firstMom_thirdFix = np.zeros(shape=(rank)), np.zeros(shape=(rank)), np.zeros(shape=(rank))
            for ins_out_pair in exampleIDs[rule]:
                in_idx = [int(num) for num in re.findall("In:(\d+)", ins_out_pair)]
                out_idx = int(re.findall("Out:(\d+)", ins_out_pair)[0])
                outerProd += tensorProduct(Z[out_idx,:], Y[in_idx[0],:], Y[in_idx[1],:])                            
                secMom_firstFree += np.outer(Y[in_idx[0], :], Y[in_idx[1], :])
                secMom_secFree += np.outer(Z[out_idx,:], Y[in_idx[1], :])
                secMom_thirdFree += np.outer(Z[out_idx,:], Y[in_idx[0], :])
                firstMom_firstFix += Z[out_idx,:]
                firstMom_secFix += Y[in_idx[0], :]
                firstMom_thirdFix += Y[in_idx[1], :]
            outerProd = np.multiply(outerProd, 1.0/ruleCount)
            firstRes = matrixVectorOuterProd(np.multiply(secMom_thirdFree, 1.0/ruleCount), np.multiply(firstMom_thirdFix, 1.0/ruleCount))
            secondRes = matrixVectorOuterProd(np.multiply(secMom_secFree, 1.0/ruleCount), np.multiply(firstMom_secFix, 1.0/ruleCount))
            thirdRes = matrixVectorOuterProd(np.multiply(secMom_firstFree, 1.0/ruleCount), np.multiply(firstMom_firstFix, 1.0/ruleCount))
            E_2 = np.multiply(firstRes, 1.0/3) + np.multiply(secondRes, 1.0/3) + np.multiply(thirdRes, 1.0/3)
            E_3 = tensorProduct(np.multiply(firstMom_firstFix, 1.0/ruleCount), np.multiply(firstMom_secFix, 1.0/ruleCount), np.multiply(firstMom_thirdFix, 1.0/ruleCount))
            K = coeff*E_3 + (1-coeff)*E_4_tensor
            param = coeff*outerProd + (1-coeff)*(coeff*E_2 + (1-coeff)*K)
        elif arity == 2:
            outerProd = np.zeros(shape=(rank, rank))
            firstMom_firstFix, firstMom_secFix = np.zeros(shape=(rank)), np.zeros(shape=(rank))
            for ins_out_pair in exampleIDs[rule]:
                in_idx = int(re.findall("In:(\d+)", ins_out_pair)[0])
                out_idx = int(re.findall("Out:(\d+)", ins_out_pair)[0])
                outerProd += np.outer(Z[out_idx,:], Y[in_idx,:])
                firstMom_firstFix += Z[out_idx, :]
                firstMom_secFix += Y[in_idx, :]
            outerProd = np.multiply(outerProd, 1.0/ruleCount)
            E_3 = np.outer(np.multiply(firstMom_firstFix, 1.0/ruleCount), np.multiply(firstMom_secFix, 1.0/ruleCount))
            K = coeff*E_3 + (1-coeff)*E_4_matrix
            param = coeff*outerProd + (1-coeff)*K
        elif arity == 1: #no smoothing for lexical rules
            outerProd = np.zeros(shape=(rank))
            for ins_out_pair in exampleIDs[rule]:                
                out_idx = int(re.findall("Out:(\d+)", ins_out_pair)[0])
                outerProd += Z[out_idx,:]
            param = np.multiply(outerProd, 1.0/ruleCount)
        else:
            sys.stderr.write("Rule has more than 2 inside trees or more than 1 outside tree!\n%s\n"%rule)
        param = np.multiply(param, scale) #scale by MLE 
        elements = rule.split(' ||| ')
        src_key = ' ||| '.join(elements[:-1])
        tgt_key = elements[-1]
        srcDict = paramDict[src_key] if src_key in paramDict else {}
        srcCountDict = countDict[src_key] if src_key in countDict else {}
        srcDict[tgt_key] = param
        srcCountDict[tgt_key] = len(exampleIDs[rule]) #for MLE estimates
        paramDict[src_key] = srcDict
        countDict[src_key] = srcCountDict
    return (paramDict, countDict)

'''
function to compute outer product between a matrix and a vector. 
'''        
def matrixVectorOuterProd(matrix, vector):
    rank = vector.shape[0]
    result = np.zeros(shape=(rank, rank, rank))
    idx = 0
    for val in np.nditer(vector):
        result[:,:,idx] = np.multiply(matrix, val)
        idx += 1
    return result

'''
the main function that computes the parameters associated with 
the various rules by computing outer products of feature vectors. 
Row indices to the inside and outside feature vectors are stored
and used to compute the outer products. We also maintain MLE counts
for each source phrase rule. 
'''
def computeCorrelations(Y, Z):    
    paramDict = {}
    countDict = {}
    numExamples = Y.shape[0]
    rank = Y.shape[1]
    scale = 1.0 / numExamples
    for rule in exampleIDs: #rule is an actual rule LHS ||| src ||| tgt
        arity = len(exampleIDs[rule][0].split())
        #scale = 1.0 / len(exampleIDs[rule])
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
        srcCountDict = countDict[src_key] if src_key in countDict else {}
        srcDict[tgt_key] = outerProd
        srcCountDict[tgt_key] = len(exampleIDs[rule])
        paramDict[src_key] = srcDict
        countDict[src_key] = srcCountDict
    return (paramDict, countDict)

'''
Function used in MLE for normalization purposes. 
'''
def normalizeCountDict(countDict):
    normalizer = 0
    for srcKey in countDict:
        #normalizer = sum([countDict[srcKey][tgtKey] for tgtKey in countDict[srcKey]])
        normalizer += sum([countDict[srcKey][tgtKey] for tgtKey in countDict[srcKey]]) #joint normalization
    for srcKey in countDict: 
        for tgtKey in countDict[srcKey]:
            countDict[srcKey][tgtKey] /= float(normalizer)
        print "Source key: %s"%(srcKey)
        print "Target distribution: "
        print countDict[srcKey]
    return countDict

'''
function to compute generalized outer product (tensor product)
for 3 vectors.  
'''        
def tensorProduct(vec1, vec2, vec3):
    vs = [vec1, vec2, vec3]
    shape = map(len, vs)
    #specify the orientation of each vector before multiplication
    newshapes = np.diag(np.array(shape)-1)+1
    reshaped = [x.reshape(y) for x,y in zip(vs, newshapes)]
    #then, multiply directly
    return reduce(np.multiply, np.ix_(*vs))

'''
given a list of dictionaries (each dictionary corresponding
to a rule, where the key is a featureID and the value is 1)
this function converts the dictionaries into a SciPy sparse
feature matrix. 
'''
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
    #scaledFeat = centerFeatures(sparseFeat)
    #scaledFeat = sparseFeat    
    return scaledFeat

def centerFeatures(sparseFeat):
    meanShift = sparseFeat - (sparseFeat.sum(axis=0) * (1.0 / sparseFeat.shape[0]))
    return meanShift

'''
this function rescales the features based on the variance of 
the feature counts, so more common features get scaled down
and less frequent features get scaled up.  
'''
def rescaleFeatures(sparseFeat, kappa):
    spFeatTransSquared = sparseFeat.transpose(copy=True)
    spFeatTransSquared.data **= 2
    scaleDenom = spFeatTransSquared.sum(axis=1) + kappa
    scaleVec = np.sqrt((sparseFeat.shape[0] - 1) * np.reciprocal(scaleDenom))
    scaleVecSp = sp.spdiags(scaleVec.flatten(), [0], len(scaleVec), len(scaleVec))
    return sparseFeat * scaleVecSp

'''
this function goes through the singleton rules (rules that
only appear once in the corpus) and eliminates them from 
the feature matrix.  Since we only look at singleton pre-terminals,
we do not need to change the features of the outside matrix. 
Once features are updated, we also update the example IDs correspondingly. 

EDIT: make changes
'''
def computeOOVProbMass(inFeatures):
    global exampleIDs
    preTermSingletons = [k for k,v in exampleIDs.items() if len(v) == 1 and len(v[0].split()) == 1] #val is a string that provides current, inside and outside tree IDs    
    sys.stderr.write("Number of pre-term singletons: %d\n"%(len(preTermSingletons)))
    sys.stderr.write("Out of %d rule types and %d rule tokens\n"%(len(exampleIDs), len(inFeatures)))
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
    sys.stderr.write("Number of pre-term singletons with valid (unary or binary) parents: %d\n"%numValidParents)
    return preTermSingletons

'''
the main feature extraction function.  A recursive function that, 
given a tree, controls the calls to the specific feature extraction
functions, and adds the results to the list of dictionaries.  It then
recursively calls update_features on the children, and also maintains
the example list (for a given unique rule, where it occurs, who its
inside and outside trees where, etc.)
'''
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
        if "lengthF" in featBinDict:
            extractSpanLengthFeatures(sent_tree, insFeatDict, outFeatDict)
        #define additional feature handling here
        inFeat.append(insFeatDict)
        outFeat.append(outFeatDict)
        curID = exampleID
        exampleID += 1
    else:
        sys.stderr.write("Rule %s has more than 2 NTs, filtered out from parameter estimation\n"%(sent_tree.rule))
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

def bucketSpanLength(span_length):
    if span_length <= 5:
        return span_length
    elif 5 < span_length and span_length <= 10:
        return 6
    elif 10 < span_length and span_length <= 20:
        return 7
    else:
        return 8

def extractSpanLengthFeatures(sent_tree, inFeatDict, outFeatDict):
    srcSpanLength = bucketSpanLength(len(sent_tree.srcYield().split()))
    tgtSpanLength = bucketSpanLength(len(sent_tree.tgtYield().split()))
    inFeatDict[addCheckFeature("srcSpanLength_%d"%srcSpanLength, "in")] = 1
    inFeatDict[addCheckFeature("tgtSpanLength_%d"%tgtSpanLength, "in")] = 1
    leftChild = True
    for child in sent_tree.children: #add length features of children
        srcSpanLengthChild = bucketSpanLength(len(child.srcYield().split()))
        tgtSpanLengthChild = bucketSpanLength(len(child.tgtYield().split()))
        if leftChild:
            inFeatDict[addCheckFeature("srcSpanLengthLeftChild_%d"%srcSpanLengthChild, "in")] = 1
            inFeatDict[addCheckFeature("tgtSpanLengthLeftChild_%d"%tgtSpanLengthChild, "in")] = 1
        else:
            inFeatDict[addCheckFeature("srcSpanLengthRightChild_%d"%srcSpanLengthChild, "in")] = 1
            inFeatDict[addCheckFeature("tgtSpanLengthRightChild_%d"%tgtSpanLengthChild, "in")] = 1
        leftChild = False
    if sent_tree.parent is not None:
        srcSpanLengthParent = bucketSpanLength(len(sent_tree.parent.srcYield().split()))
        tgtSpanLengthParent = bucketSpanLength(len(sent_tree.parent.tgtYield().split()))
        outFeatDict[addCheckFeature("srcSpanLengthParent_%d"%srcSpanLengthParent, "out")] = 1
        outFeatDict[addCheckFeature("tgtSpanLengthParent_%d"%tgtSpanLengthParent, "out")] = 1
        leftChild = True
        for child in sent_tree.parent.children:
            if child is not sent_tree:
                srcSpanLengthSibling = bucketSpanLength(len(child.srcYield().split()))
                tgtSpanLengthSibling = bucketSpanLength(len(child.tgtYield().split()))
                if leftChild:
                    outFeatDict[addCheckFeature("srcSpanLengthLeftSibling_%d"%srcSpanLengthSibling, "out")] = 1
                    outFeatDict[addCheckFeature("tgtSpanLengthLeftSibling_%d"%tgtSpanLengthSibling, "out")] = 1
                else:
                    outFeatDict[addCheckFeature("srcSpanLengthLeftSibling_%d"%srcSpanLengthSibling, "out")] = 1
                    outFeatDict[addCheckFeature("tgtSpanLengthLeftSibling_%d"%tgtSpanLengthSibling, "out")] = 1
                leftChild = False

'''
extracts words (or word classes) of any lexical items
in the current rule, as well as those in the yield of the
left and right children, and the outside tree. We currently
look at only the first and last words in the yield. 
'''
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
