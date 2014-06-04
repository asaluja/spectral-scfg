#!/usr/bin/python -tt

'''
File: em_estimation.py
Author: Avneesh Saluja (avneesh@cs.cmu.edu)
Date: May 11, 2014
Description: main script for EM estimation of parameters. 
Write more detailed explanation here. 
'''

import sys, commands, string, gzip, getopt, cPickle, os, re, math, Queue, time
import numpy as np
#np.seterr(under='raise')
import multiprocessing as mp
from tree import tree

params = {} #current set of parameters is global
counts = {}
epsilon = 1e-15
scaling = 1
'''
Given 3 vectors, computes tensor product (generalized 3-way outer product)
'''
def TensorProduct(vec1, vec2, vec3):
    vs = [vec1, vec2, vec3]
    shape = map(len, vs)
    #specify the orientation of each vector before multiplication
    newshapes = np.diag(np.array(shape)-1)+1
    reshaped = [x.reshape(y) for x,y in zip(vs, newshapes)]
    #then, multiply directly
    return reduce(np.multiply, np.ix_(*vs))

'''
simple function to increment count of rule in counts dictionary.
This dictionary tells us which rules are singletons (for OOV filtering). 
'''
def AddCounts(rule):
    curCount = counts[rule] if rule in counts else 0
    counts[rule] = curCount + 1

'''
this function normalizes the parameter structure. 
If it is a binary rule, then it sums over h2 and h3 (the 2nd and 3rd
modes of the tensor) and normalizes for every h1. 
Similarly, for unary rules with matrices, we sum over all columns
for a given row and normalize, and for vectors we just normalize by
the sum of the vector. 
NOTE: we are not using this function anymore, deprecated
'''
def NormalizeParam(unnormalized):
    order = len(unnormalized.shape)
    rank = unnormalized.shape[0]    
    normalized = None
    if order == 3:
        normalized = np.zeros(shape=(rank, rank, rank))
        for h1 in xrange(0, rank):
            normalizer = np.sum(unnormalized[h1, :, :])
            normalized[h1, :, :] = np.divide(unnormalized[h1, :, :], normalizer)
    elif order == 2:
        normalized = np.zeros(shape=(rank, rank))
        for h1 in xrange(0, rank):
            normalizer = np.sum(unnormalized[h1, :])
            normalized[h1, :] = np.divide(unnormalized[h1, :], normalizer)
    elif order == 1:
        normalized = np.zeros(shape=(rank))
        normalizer = np.sum(unnormalized)
        normalized = np.divide(unnormalized, normalizer)
    return normalized

'''
Simple function that returns the arity of a rule. 
'''
def CheckArity(rule):
    srcRule = rule.split(' ||| ')[1]
    return len(re.findall(r'\[([^]]+)\]', srcRule))

'''
Simple function that does a top-down traversal of a 
training tree and checks if all rules are valid (arity <= 2)
'''
def CheckForInvalidRules(training_tree):
    if len(training_tree.children) <= 2: 
        valid = True
        for child in training_tree.children:
            valid = valid & CheckForInvalidRules(child)
        return valid
    else:
        return False

'''
this function is used during initialization; we go through all of
the parameters in the training corpus, and initialize the appropriate
structure (tensor for binary rules, matrix for unary rules, and vector
for preterminal rules) and make sure the structure is appropriately normalized. 
'''
def InitCounts(training_tree, rank):
    AddCounts(training_tree.rule)
    if training_tree.rule not in params:        
        arity = CheckArity(training_tree.rule)        
        unnormalized = None
        if arity == 2:
            unnormalized = np.random.random_sample((rank, rank, rank))
        elif arity == 1:
            unnormalized = np.random.random_sample((rank, rank))
        elif arity == 0:
            unnormalized = np.random.random_sample((rank,))
        else: #filtered out invalid trees when reading in, so this should never fire 
            sys.stderr.write("Rule %s has more than 2 NTs, not included in parameter dictionary\n"%training_tree.rule)
        if unnormalized is not None:
            params[training_tree.rule] = unnormalized 
            #params[training_tree.rule] = NormalizeParam(unnormalized)
    for child in training_tree.children: #at this stage, recurse on children
        InitCounts(child, rank)

'''
this function implements Matsuzaki initialization.
Explain Matsuzaki initialization here. 
'''
def InitCountsMatsuzaki(training_tree, rank, mle_params):
    a = -math.log(3)
    b = math.log(3)
    AddCounts(training_tree.rule)
    if training_tree.rule not in params:
        arity = CheckArity(training_tree.rule)
        rule = training_tree.rule.split(' ||| ')
        srcKey = ' ||| '.join(rule[0:2])
        tgtKey = rule[2]
        if srcKey in mle_params and tgtKey in mle_params[srcKey]: #only add to params if also in MLE params
            mle = mle_params[srcKey][tgtKey]
            unnormalized = None
            if arity == 2:
                unnormalized = (b-a)*np.random.random_sample((rank, rank, rank))+a
            elif arity == 1:
                unnormalized = (b-a)*np.random.random_sample((rank, rank))+a
            elif arity == 0:
                unnormalized = (b-a)*np.random.random_sample((rank,))+a
            else:
                sys.stderr.write("Rule %s has more than 2 NTs, not included in parameter dictionary\n"%training_tree.rule)
            if unnormalized is not None:
                unnormalized = np.exp(unnormalized)
                unnormalized *= mle
                params[training_tree.rule] = unnormalized
                #params[training_tree.rule] = NormalizeParam(unnormalized)
    for child in training_tree.children:
        InitCountsMatsuzaki(child, rank, mle_params)

def ComputeInside(training_tree, OOV, alpha): 
    alpha_children = []
    for child in training_tree.children:
        alpha_children.append(ComputeInside(child, OOV, alpha))
    arity = len(training_tree.children)
    if arity == 0:
        alpha[training_tree] = params["[X] ||| <unk> ||| <unk>"] if OOV and counts[training_tree.rule] == 1 else params[training_tree.rule]
    elif arity == 1:
        alpha[training_tree] = params[training_tree.rule].dot(alpha_children[0]); 
    elif arity == 2:
        x1_alpha = alpha_children[0]
        x2_alpha = alpha_children[1]
        result = np.tensordot(x2_alpha, params[training_tree.rule], axes=[0,2])
        alpha[training_tree] = result.dot(x1_alpha)
    return alpha[training_tree]

def ComputeOutside(training_tree, alpha, beta):
    arity = len(training_tree.children)    
    left = True
    for child in training_tree.children:
        if arity == 1:
            beta[child] = beta[training_tree].dot(params[training_tree.rule])
        elif arity == 2:
            result = np.tensordot(beta[training_tree], params[training_tree.rule], axes=[0,0])
            if left:
                beta[child] = result.dot(alpha[training_tree.children[1]])
                left = False
            else:
                beta[child] = alpha[training_tree.children[0]].dot(result)
    for child in training_tree.children:
        ComputeOutside(child, alpha, beta)

def ConvertToTripleDict(training_tree, OOV, alpha, beta, ab_dict):
    arity = len(training_tree.children)
    vectors = None
    if arity == 0:
        vectors = (beta[training_tree],)
    elif arity == 1:
        vectors = (beta[training_tree], alpha[training_tree.children[0]])
    elif arity == 2:
        vectors = (beta[training_tree], alpha[training_tree.children[0]], alpha[training_tree.children[1]])
    if training_tree.parent is None: #i.e., root level rule
        ab_dict["Root"] = alpha[training_tree] #contains inside probability of entire tree
    rule = "[X] ||| <unk> ||| <unk>" if OOV and arity == 0 and counts[training_tree.rule] == 1 else training_tree.rule
    if rule in ab_dict:
        ab_dict[rule].append(vectors)
    else:
        ab_dict[rule] = [vectors]
    for child in training_tree.children:
        ConvertToTripleDict(child, OOV, alpha, beta, ab_dict)

def InsideOutside(training_examples, OOV, out_q):
    ins_out_probs = []
    for training_tree in training_examples:
        alpha = {}
        ComputeInside(training_tree, OOV, alpha) 
        beta = {}
        beta[training_tree] = params["Pi"]
        ComputeOutside(training_tree, alpha, beta)
        alpha_beta_dict = {}
        ConvertToTripleDict(training_tree, OOV, alpha, beta, alpha_beta_dict)
        ins_out_probs.append(alpha_beta_dict)
    out_q.put(ins_out_probs)

'''
The inside/outside probability vectors for each non-terminal node for each sentence's
minimal tree is stored in ab_list (indexed by sentence, each entry is a dict).  
'''
def UpdateParameters(ab_list):
    new_params = {}
    ll = 0
    count = 0
    for tree_dict in ab_list: #indexed per sentence
        g = tree_dict["Root"].dot(params["Pi"])        
        ll += math.log(g)
        count += 1
        if g == 0:
            print "normalizer = 0 for example number %d"%count
        for rule in tree_dict: #for each minimal rule for the sentence, have inside/outside probs of spans
            if rule != "Root":
                arity = CheckArity(rule)
                for vectors in tree_dict[rule]:
                    result = None
                    if arity == 0:
                        result = vectors[0]
                    elif arity == 1:
                        result = np.outer(vectors[0], vectors[1])
                    elif arity == 2:
                        result = TensorProduct(vectors[0], vectors[1], vectors[2])
                    else:
                        sys.stderr.write("Fatal error! Found rule with arity > 2 while updating parameters\n")
                        sys.exit()
                    result = np.multiply(params[rule], result)
                    result = np.divide(result, g)
                    if rule in new_params:
                        new_params[rule] += result
                    else:
                        new_params[rule] = result
    print "Log Likelihood: %.3f"%ll
    return new_params

def UpdatePiParameters(ins_out_probs, new_params):
    param = params["Pi"]
    new_pi = np.zeros(param.shape)
    for alpha_beta_dict in ins_out_probs:
        alpha_top = alpha_beta_dict["Root"]
        g = alpha_top.dot(param)
        new_pi += np.divide(np.multiply(alpha_top, param), g)
    new_params["Pi"] = new_pi

def ExpectationStep(training_examples, OOV, numProc):
    out_q = mp.Queue()
    chunksize = int(math.ceil(len(training_examples) / float(numProc)))
    #InsideOutside(training_examples, out_q)
    #ins_out_probs = out_q.get()
    procs = []
    for i in range(numProc):
        end = len(training_examples) if i == numProc - 1 else (i+1)*chunksize
        p = mp.Process(target=InsideOutside,args=(training_examples[chunksize*i : end], OOV, out_q))
        procs.append(p)
        p.start()
    ins_out_probs = []
    for i in range(numProc): #collect all results in a single list
        ins_out_probs += out_q.get()
    for p in procs:
        p.join()        
    new_params = UpdateParameters(ins_out_probs)
    UpdatePiParameters(ins_out_probs, new_params)
    return new_params

def MaximizationStep(paramDict, rank):
    normalizer = np.zeros((rank,))
    normalizer_pi = np.zeros((rank,))
    for rule in paramDict: #first loop to compute the normalizer for a fixed h1
        order = len(paramDict[rule].shape)
        if rule == "Pi": #normalize separately
            normalizer_pi = np.sum(paramDict[rule])
        else:
            for h1 in xrange(0, rank):
                if order == 3:
                    normalizer[h1] += np.sum(paramDict[rule][h1, :, :])
                elif order == 2:
                    normalizer[h1] += np.sum(paramDict[rule][h1,:])
                elif order == 1:
                    normalizer[h1] += paramDict[rule][h1]
    for rule in paramDict: #then loop through and rescale each parameter
        order = len(paramDict[rule].shape)
        if rule == "Pi":
            paramDict[rule] = np.divide(paramDict[rule], normalizer_pi)
        else:
            for h1 in xrange(0, rank):
                if order == 3:
                    paramDict[rule][h1,:,:] = np.divide(paramDict[rule][h1,:,:], normalizer[h1])                
                elif order == 2:
                    paramDict[rule][h1,:] = np.divide(paramDict[rule][h1,:], normalizer[h1])
                elif order == 1:
                    paramDict[rule][h1] = np.divide(paramDict[rule][h1], normalizer[h1])
        paramDict[rule] = np.multiply(paramDict[rule], scaling)

def AddOOVParameter(rank, matsu, mle_params):
    if matsu:
        srcKey = "[X] ||| <unk>"
        tgtKey = "<unk>"
        mle_param = mle_params[srcKey][tgtKey]
        a = -math.log(3)
        b = math.log(3)    
        unnormalized = (b-a)*np.random.random_sample((rank,))+a
        unnormalized = np.exp(unnormalized)
        unnormalized *= mle_param
    else:
        unnormalized = np.random.random_sample((rank,))
    params["[X] ||| <unk> ||| <unk>"] = unnormalized
    #params["[X] ||| <unk> ||| <unk>"] = NormalizeParam(unnormalized)

def WriteOutParameters(filename, isFilter, filteredRules):
    param_writeout = {}
    num_params = 0
    preterm_singletons = [k for k,v in counts.items() if v == 1 and CheckArity(k) == 0] 
    OOV_param = params["[X] ||| <unk> ||| <unk>"]
    avgOOV_param = np.divide(OOV_param, len(preterm_singletons)+1)
    avgOOV_param = np.divide(avgOOV_param, scaling)
    for param in params: #unk not in params??
        temp_param = np.divide(params[param], scaling)
        if param != "Pi":
            srcKey = ' ||| '.join(param.split(' ||| ')[:-1])
            tgtKey = param.split(' ||| ')[-1]
            srcDict = param_writeout[srcKey] if srcKey in param_writeout else {}
            if isFilter:
                if srcKey in filteredRules and tgtKey in filteredRules[srcKey]:
                    srcDict[tgtKey] = avgOOV_param if param == "[X] ||| <unk> ||| <unk>" or param in preterm_singletons else temp_param
                    num_params += 1
                    param_writeout[srcKey] = srcDict
            else:
                srcDict[tgtKey] = avgOOV_param if param == "[X] ||| <unk> ||| <unk>" or param in preterm_singletons else temp_param
                num_params += 1
                param_writeout[srcKey] = srcDict
        else:
            param_writeout["Pi"] = params[param]
            num_params += 1
    print "Writing out %d parameters (after potential filtering)"%num_params
    cPickle.dump(param_writeout, open(filename, "wb"))
            
def main():
    global params, counts, scaling
    optsDict = {}
    (opts, args) = getopt.getopt(sys.argv[1:], 'f:m:n:o')
    for opt in opts:
        if opt[0] == '-f': #filter rules
            optsDict["filter"] = opt[1]
        elif opt[0] == '-n': #number of processes
            optsDict["numProc"] = int(opt[1])
        elif opt[0] == '-m': #matsuzaki initialization
            optsDict["Matsuzaki"] = opt[1] #contains location of MLE params
        elif opt[0] == '-o': #OOV
            optsDict["OOV"] = 1
    minrule_grammars_loc = args[0]
    rank = int(args[1])
    scaling = int(args[4])
    if "numProc" not in optsDict:
        optsDict["numProc"] = 4 #default value for numProc

    numSentences = len(os.listdir(minrule_grammars_loc))
    training_examples = [] #list of trees
    start = time.clock()
    for lineNum in range(0, numSentences): #read in the training trees in the beginning into a dict
        minrule_fh = gzip.open(minrule_grammars_loc + "grammar.%d.gz"%(lineNum))
        sync_tree = tree(0, None, None, minrule_fh)
        if CheckForInvalidRules(sync_tree):
            training_examples.append(sync_tree)
        else:
            sys.stderr.write("Training example %d has invalid tree (contains rule with # NTs on RHS > 2)\n"%lineNum)
    timeTaken = time.clock() - start
    print "Read in %d training tree examples.  Time taken: %.3f sec"%(len(training_examples), timeTaken)    

    params = {}
    mle_params = None
    #np.random.seed(42)
    if "Matsuzaki" in optsDict:
        params_fh = open(optsDict["Matsuzaki"], 'rb')
        mle_params = cPickle.load(params_fh)
    start = time.clock()
    for training_tree in training_examples: #only loops through valid rules
        if "Matsuzaki" in optsDict:
            InitCountsMatsuzaki(training_tree, rank, mle_params)
        else:
            InitCounts(training_tree, rank)
    #params["Pi"] = NormalizeParam(np.random.random_sample((rank,)))    
    params["Pi"] = np.random.random_sample((rank,))
    if "OOV" in optsDict:
        AddOOVParameter(rank, "Matsuzaki" in optsDict, mle_params)
        for param in params.keys():
            if param != "Pi" and param != "[X] ||| <unk> ||| <unk>" and CheckArity(param) == 0 and counts[param] == 1:
                del params[param]
    else:
        params["[X] ||| <unk> ||| <unk>"] = np.zeros((rank,))
    MaximizationStep(params, rank) #normalize params here
    timeTaken = time.clock() - start    
    print "Initialized %d parameters for rules in training corpus.  Time taken: %.3f sec"%(len(params), timeTaken)

    filteredRules = None
    if "filter" in optsDict:
        filteredRules = cPickle.load(open(optsDict["filter"], 'rb'))
    numIters = int(args[2])
    outDir = args[3]
    start = time.clock()
    start_iter = start
    for iterNum in range(0, numIters):        
        new_params = ExpectationStep(training_examples, "OOV" in optsDict, optsDict["numProc"])
        print "Iteration %d: E-step completed; computed partial counts and updated parameters"%(iterNum+1)
        MaximizationStep(new_params, rank)
        print "Iteration %d: M-step completed; normalized parameters"%(iterNum+1)
        params = new_params
        if "OOV" not in optsDict:
            params["[X] ||| <unk> ||| <unk>"] = np.zeros((rank,))
        print "Number of params in new params: %d"%len(new_params)
        if (iterNum + 1)%5 == 0: #write out every 10 iterations
            WriteOutParameters(outDir + "/em.iter%d.params"%(iterNum+1), "filter" in optsDict, filteredRules)
        timeTaken = time.clock() - start_iter
        print "Iteration %d complete. Time taken: %.3f sec"%(iterNum+1, timeTaken)
        start_iter = time.clock()
    print "EM complete.  Time taken: %.3f sec"%(time.clock()-start)
    
    
if __name__ == "__main__":
    main()
