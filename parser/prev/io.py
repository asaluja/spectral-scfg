#!/usr/bin/python -tt

'''
File: inside_outside.py
Date: November 17, 2013
Description: this script takes as arguments
a cPickle file with the latent-SCFG parameters
in a dictionary (key: the rule ID) and a corpus
in stdin, and computes marginals \mu(X,i,j) for
1 \leq i \leq j \eq N, where N is the length of
the sentence.  These marginals are printed to stdout.
Arguments:
arg1: dictionary of parameters; key is a src-tgt phrase pair in tuple format; value is parameters (tensor, matrix, or vector)
arg2: rank of latent feature space
arg3: maximum number of lexical items in a rule
'''

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
    if arity in vals:
        vals[arity].extend(cell_rules)
    else:
        vals[arity] = cell_rules
    rules_used[key] = vals    

def matchPreTerminals(words):    
    paramDictKeys = [pp for pp in paramDict if pp != "Pi"]
    pps = [pp for pp in paramDictKeys if pp.split(' ||| ')[1] == words]
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
        paramDictKeys = [pp for pp in paramDict if pp != "Pi"]
        pps.extend([(pp, lex_lens[0], lex_lens[1]) for pp in paramDictKeys if pp.split(' ||| ')[1] == srcPhr.lstrip().rstrip()])
    pps = [pp for pp in pps if len(pp) > 0] #filter empty entries out
    return pps

def computeUnaryRulesAlpha(words, start, end, alpha):
    pps = matchUnaryRules(words, start, end)    
    aggregate = np.zeros((rank))
    rules = []
    for pp in pps:
        (key, f1_len, f2_len) = pp
        aggregate += paramDict[key].dot(alpha[start+f1_len, end-start-f2_len-f1_len-1, :])
        rules.append((key, (start+f1_len, end-start-f2_len-f1_len-1)))
    addPPsToRulesUsed(start, end, "arity1", rules)
    return aggregate

def computeUnaryRulesBeta(cell_rules, beta, curCell_idx):
    for rule in cell_rules:
        key, cell_idx = rule
        beta[cell_idx[0], cell_idx[1]] += beta[curCell_idx[0], curCell_idx[1]].dot(paramDict[key])

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
        paramDictKeys = [pp for pp in paramDict if pp != "Pi"]
        pps.extend([(pp, lex_lens[0], lex_lens[1], lex_lens[2]) for pp in paramDictKeys if pp.split(' ||| ')[1] == srcPhr.lstrip().rstrip()])
    pps = [pp for pp in pps if len(pp) > 0]
    return pps

def computeBinaryRulesAlpha(words, start, end, split, alpha):
    pps = matchBinaryRules(words, start, end, split)    
    aggregate = np.zeros((rank))
    rules = []
    for pp in pps:
        (key, f1_len, f2_len, f3_len) = pp        
        x1_alpha = alpha[start+f1_len, split-start-f1_len, :]
        x2_alpha = alpha[split+1+f2_len, end-split-2-f2_len-f3_len, :]
        result = np.tensordot(paramDict[key], x1_alpha, axes=([1], [0]))
        result = np.tensordot(result, x2_alpha, axes=([1], [0]))
        aggregate += result
        rules.append((key, (start+f1_len, split-start-f1_len), (split+1+f2_len, end-split-2-f2_len-f3_len)))
    addPPsToRulesUsed(start, end, "arity2", rules)
    return aggregate

def computeBinaryRulesBeta(cell_rules, alpha, beta, curCell_idx):
    for rule in cell_rules:        
        key, cell_idx_left, cell_idx_right = rule
        x_beta = beta[curCell_idx[0], curCell_idx[1]]
        result = np.tensordot(paramDict[key], x_beta, axes=([0], [0]))
        x_alpha_right = alpha[cell_idx_right[0], cell_idx_right[1]]
        x_alpha_left = alpha[cell_idx_left[0], cell_idx_left[1]]
        right_beta = np.tensordot(result, x_alpha_left, axes=([0], [0]))
        beta[cell_idx_right[0], cell_idx_right[1]] += right_beta
        left_beta = np.tensordot(result, x_alpha_right, axes=([1], [0]))
        beta[cell_idx_left[0], cell_idx_left[1]] += left_beta

def computeAlpha(words):
    N = len(words)
    alpha = np.zeros((N, N, rank))
    for i in range(0, N): #base case initialization
        alpha[i, 0, :] += computePreTerms(words, i, i+1)
    for j in range(2, N+1): #length of span from 2 to N 
        for i in range(0, N-j+1): #left (start) of span            
            if j <= longest_lex: #only check for pre-term parameters if length of pre-term does not exceed max pre-term length
                alpha[i,j-1,:] += computePreTerms(words, i, i+j)            
            alpha[i,j-1,:] += computeUnaryRulesAlpha(words, i, i+j, alpha) #don't need split point for unary rules       
            for k in range(i, i+j-1): 
               alpha[i,j-1,:] += computeBinaryRulesAlpha(words, i, i+j, k, alpha)
    return alpha

def computeBeta(N, alpha):
    beta = np.zeros((N, N, rank)) 
    beta[0,N-1,:] = paramDict['Pi'] #base case
    for j in range(N, 1, -1): #length of span from entire sentence to just 1 word
        for i in range(0, N-j+1): #start (left) of span
            cell_rules = rules_used[(i,j-1)]                           
            if "arity1" in cell_rules:
                computeUnaryRulesBeta(cell_rules["arity1"], beta, (i,j-1))
            if "arity2" in cell_rules:
                computeBinaryRulesBeta(cell_rules["arity2"], alpha, beta, (i, j-1))
    return beta

def computeMarginals(words):
    alpha = computeAlpha(words) #N x N x m tensor, where N is sentence length and m is rank        
    beta = computeBeta(len(words), alpha)
    marginals = np.zeros((len(words), len(words)))
    for i in range(0, len(words)):
        for j in range(0, len(words)):
            marginals[i,j] = np.dot(alpha[i,j,:], beta[i,j,:])
            print "cell item [%d,%d] alpha, beta, and marginal"%(i,j)
            print alpha[i,j,:]
            print beta[i,j,:]
            print marginals[i,j]
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
        #for i in range(0,N):
        #    for j in range(0,N):
        #        print marginals[i,j],
        #    print
        sys.stdout.flush()
        sentence_marginals.append(marginals)
    #write out marginals to file here

if __name__ == "__main__":
    main()
