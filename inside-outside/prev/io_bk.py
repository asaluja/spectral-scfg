#!/usr/bin/python -tt

###########################
#File: inside_outside.py
#Date: November 17, 2013
#Description: this script takes as arguments
#a cPickle file with the latent-SCFG parameters
#in a dictionary (key: the rule ID) and a corpus
#in stdin, and computes marginals \mu(X,i,j) for
#1 \leq i \leq j \eq N, where N is the length of
#the sentence.  These marginals are printed to stdout.
#Arguments:
#arg1: dictionary of parameters; key is a src-tgt phrase pair in tuple format; value is parameters (tensor, matrix, or vector)
#arg2: rank of latent feature space
#arg3: maximum number of lexical items in a rule
###########################

import sys, commands, string, cPickle, itertools
import numpy as np

params_fh = open(sys.argv[1], 'rb')
paramDict = cPickle.load(params_fh)
rank = int(sys.argv[2])
longest_lex = int(sys.argv[3])
rules_used = {} #global data structure that stores rules used during inside computation, to make outside computation easier

def addPPsToRulesUsed(start, end, arity, cell_rules):
    global rules_used
    key = (start, end-start-1) #store key in parent
    vals = rules_used[key] if key in rules_used else {}
    vals[arity] = cell_rules
    rules_used[key] = vals    

def matchPreTerminals(words):    
    pps = [pp for pp in paramDict if pp[0] == words]
    return pps

def computePreTerms(words, start, end):
    span = ' '.join(words[start:end])    
    pps = matchPreTerminals(span)
    aggregate = np.zeros((rank))
    for pp in pps:
        aggregate += paramDict[pp] 
    return aggregate

def sum_to_n(n, size):
    if size == 1:
        yield [n]
        return
    start = (n + size - 1) // size
    stop = n
    for i in range(start, stop + 1):
        for tail in sum_to_n(n - i, size - 1):
            yield [i] + tail

def genLexLengths(span, lengths, arity):
    if span == 0:
        return
    else:
        for partition in sum_to_n(span, arity):
            lengths.extend(list(itertools.permutations(partition)))
        genLexLengths(span-1, lengths, arity)

def matchUnaryRules(words, start, end):
    lex_len_combos = []
    span = longest_lex if (end - start - 1) > longest_lex else end - start - 1 #takes into account pre-existing info to make code run faster; -1 because we start with the case where NTs span 1 
    genLexLengths(span, lex_len_combos, 2)
    lex_len_combos = list(set(lex_len_combos))
    pps = []    
    for lex_lens in lex_len_combos:
        f1 = ' '.join(words[start:(start+lex_lens[0])])
        f2 = ' '.join(words[(end-lex_lens[1]):end])
        srcPhr = f1 + " [X,1] " + f2 #might want to change this as well
        pps.extend([(pp, lex_lens[0], lex_lens[1]) for pp in paramDict if pp[0] == srcPhr.lstrip().rstrip()])
    pps = [pp for pp in pps if len(pp) > 0] #filter empty entries out
    return pps

def computeUnaryRulesAlpha(words, start, end, alpha):
    pps = matchUnaryRules(words, start, end)    
    aggregate = np.zeros((rank))
    rules_used = []
    for pp in pps:
        (key, f1_len, f2_len) = pp
        aggregate += paramDict[key].dot(alpha[start+f1_len, end-start-f2_len-f1_len-1, :])
        rules_used.append((key, (start+f1_len, end-start-f2_len-f1_len-1)))
    addPPsToRulesUsed(start, end, "arity1", rules_used)
    return aggregate

def computeUnaryRulesBeta(words, start, end, split, beta, pps):
    aggregate = np.zeros((rank))
    for pp in pps:
        (key, f1_len, f2_len) = pp
        if split < start:
            aggregate += beta[split,end-split-1,:].dot(paramDict[key])
        else:
            aggregate += beta[start,split-start,:].dot(paramDict[key])
    return aggregate

def matchBinaryRules(words, start, end, split):
    lex_len_combos = []
    span = longest_lex if (end - start - 2) > longest_lex else end - start - 2 #takes into account pre-existing info to make code run faster; -2 because we start with case where NTs span 1 word each
    genLexLengths(span, lex_len_combos, 3)
    lex_len_combos = list(set(lex_len_combos))
    lex_len_combos = [ll for ll in lex_len_combos if ll[0] <= split - start and ll[1] + ll[2] <= end - split -2] #filtering invalid rules based on split point
    lex_len_combos = [ll for ll in lex_len_combos if ll[1] != 0] #filtering rules based on fact that 2 NTs cannot appear consecutively
    lex_len_combos.append((0,0,0))
    pps = []
    for lex_lens in lex_len_combos:
        f1 = ' '.join(words[start:(start+lex_lens[0])])
        f2 = ' '.join(words[split+1:(split+1+lex_lens[1])]) #needs to take into account arity? 
        f3 = ' '.join(words[end-lex_lens[2]+1:end])
        srcPhrList = filter(None, [f1, "[X,1]", f2, "[X,2]",f3])
        srcPhr = ' '.join(srcPhrList)
        pps.extend([(pp, lex_lens[0], lex_lens[1], lex_lens[2]) for pp in paramDict if pp[0] == srcPhr.lstrip().rstrip()])
    pps = [pp for pp in pps if len(pp) > 0]
    return pps

def computeBinaryRulesAlpha(words, start, end, split, alpha):
    pps = matchBinaryRules(words, start, end, split)    
    aggregate = np.zeros((rank))
    rules_used = []
    for pp in pps:
        (key, f1_len, f2_len, f3_len) = pp        
        x1_alpha = alpha[start+f1_len, split-start-f1_len, :]
        x2_alpha = alpha[split+1+f2_len, end-split-2-f2_len-f3_len, :]
        result = np.tensordot(paramDict[key], x1_alpha, axes=([1], [0]))
        result = np.tensordot(result, x2_alpha, axes=([1], [0]))
        aggregate += result
        rules_used.append((key, (start+f1_len, split-start-f1_len), (split+1+f2_len, end-split-2-f2_len-f3_len)))
    addPPsToRulesUsed(start, end, "arity2", rules_used)
    return aggregate

def computeBinaryRulesBeta(words, start, end, split, alpha, beta, pps):
    aggregate = np.zeros((rank))
    for pp in pps:
        (key, f1_len, f2_len, f3_len) = pp
        if split < start:
            x_beta = beta[split,end-split-1,:]
            x_alpha=alpha[split,start-split-1,:]
            result = np.tensordot(paramDict[key], x_alpha, axes=([1], [0]))
            result = np.tensordot(result, x_beta, axes=([0], [0]))
            aggregate += result        
        else:
            x_beta = beta[start,split-start,:]
            x_alpha = alpha[end,split-end,:]
            result = np.tensordot(paramDict[key], x_alpha, axes=([2], [0]))
            result = np.tensordot(result, x_beta, axes=([0], [0]))
            aggregate += result
    return aggregate

def computeAlpha(words, rules_used):
    N = len(words)
    alpha = np.zeros((N, N, rank))
    for i in range(0, N): #base case initialization
        alpha[i, 0, :] += computePreTerms(words, i, i+1, rules_used)
    for j in range(2, N+1): #length of span from 2 to N 
        for i in range(0, N-j+1): #left (start) of span            
            if j <= longest_lex: #only check for pre-term parameters if length of pre-term does not exceed max pre-term length
                alpha[i,j-1,:] += computePreTerms(words, i, i+j, rules_used)            
            alpha[i,j-1,:] += computeUnaryRulesAlpha(words, i, i+j, alpha, rules_used) #don't need split point for unary rules                        
            for k in range(i, i+j-1): 
               alpha[i,j-1,:] += computeBinaryRulesAlpha(words, i, i+j, k, alpha, rules_used)
    return alpha

def filterUnaryRules(rules, start, end, split):
    if split < start: #then current NT is right child
        f1_len = start-split #determine length allocated to outside tree outside of NT
        rules = [rule for rule in rules if rule[1] == f1_len and rule[2] == 0]
    else:
        f2_len = split-end+1
        rules = [rule for rule in rules if rule[2] == f2_len and rule[1] == 0]
    return rules

def filterBinaryRules(rules, start, end, split):
    if split < start: #then current NT is right child
        f1_f2_len = start-split-1 #-1 because the other NT has to take up a span of at least 1
        rules = [rule for rule in rules if rule[1]+rule[2] == f1_f2_len and rule[3] == 0]
    else:
        f2_f3_len = split-end
        rules = [rule for rule in rules if rule[2]+rule[3] == f2_f3_len and rule[1] == 0]
    return rules

def computeBeta(words, alpha):
    N = len(words)
    beta = np.zeros((N, N, rank)) 
    beta[0,N-1,:] = paramDict['S'] #base case
    for j in range(N-1,0,-1): #length of span from N to 2
        for i in range(0, N-j+1): #start (left) of span            
            for k in range(0, i): #acting as if it were right child
                rules = rules_used[(k,i+j-k-1)] 
                print "For tuple (i,j,k) = (%d,%d,%d), use rules from cell (%d,%d)"%(i,j,k,k,i+j-k-1)
                if "arity1" in rules:
                    filtered_rules = filterUnaryRules(rules["arity1"], i, i+j, k)
                    beta[i,j-1,:] += computeUnaryRulesBeta(words, i, i+j, k, beta, filtered_rules)
                if "arity2" in rules:
                    filtered_rules = rules["arity2"]
                    beta[i,j-1,:] += computeBinaryRulesBeta(words, i, i+j, k, alpha, beta, filtered_rules)
            for k in range(i+j, N): #acting as if it were left child
                rules = rules_used[(i,k-i)] #filter rules based on value of k
                print "For tuple (i,j,k) = (%d,%d,%d), use rules from cell (%d,%d)"%(i,j,k,i,k-i)
                if "arity1" in rules:
                    filtered_rules = filterUnaryRules(rules["arity1"], i, i+j, k)
                    beta[i,j-1,:] += computeUnaryRulesBeta(words, i, i+j, k, beta, filtered_rules)                    
                if "arity2" in rules:
                    filtered_rules = rules["arity2"]
                    beta[i,j-1,:] += computeBinaryRulesBeta(words, i, i+j, k, alpha, beta, filtered_rules)
    return beta

def computeMarginals(words):
    alpha = computeAlpha(words) #N x N x m tensor, where N is sentence length and m is rank        
    print rules_used
    beta = computeBeta(words, alpha)
    marginals = np.zeros((len(words), len(words)))
    for i in range(0, len(words)):
        for j in range(0, len(words)):
            marginals[i,j] = np.dot(alpha[i,j,:], beta[i,j,:])
            if marginals[i,j] > 1 or marginals[i,j] < 0:
                sys.stderr.write("Error! marginal of span [%d,%d] outside of range: %s\n"%(i,j,str(marginals[i,j])))
    return marginals

def main():
    sentence_marginals = [] #list of matrices
    for line in sys.stdin:
        print "Computing Chart for input sentence: "
        print line.strip()
        words = line.strip().split()
        N = len(words)
        marginals = computeMarginals(words)
        for i in range(0,N):
            for j in range(0,N):
                print marginals[i,j],
            print
        sentence_marginals.append(marginals)
    #write out marginals to file here

if __name__ == "__main__":
    main()
