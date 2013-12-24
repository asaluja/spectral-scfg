#!/usr/bin/python -tt

'''
File: hg_io.py
Date: December 20, 2013
Description: given a hypergraph representation, this file contains
the main functions to compute the inside (alpha) and outside (beta)
terms over the hypergraph.  We eventually return a dictionary to
the calling function (in intersect_scfg.py) a dictionary where the
key is the node ID, and the value is the marginal term.  This is 
then used when converting the hypergraph to a grammar. 
The current interpretation that I have is that marginals should be
defined over spans, not rules? I'm not quite sure how we can compute
marginals by multiplying alphas over tail nodes when the vector
multiplications don't work out.  
'''

import sys, commands, string, re
from trie import HyperGraph
import numpy as np

def checkArity(rule):
    return len(re.findall(r'\[([^]]+)\]', rule.split(' ||| ')[1]))

def computeInside(hg, paramDict, rank):
    PDKeys = [pp for pp in paramDict if pp != "Pi"] #filtering sentence start parameters
    alphaDict = {}
    for node in hg.nodes_: #nodes are added in topological order, so hg.nodes_ is guaranteed to be sorted
        aggregate = np.zeros((rank))
        for inEdgeID in node.in_edges_: #for each incoming hyperedge (i.e., rule that fits the span)
            src_rule = hg.edges_[inEdgeID].rule
            tail = hg.edges_[inEdgeID].tailNodes
            matching_rules = [rule for rule in PDKeys if rule.split(' ||| ')[1] == src_rule]
            for key in matching_rules:
                arity = checkArity(key)
                assert len(tail) == arity
                if arity == 0: #pre-terminal
                    aggregate += paramDict[key]
                elif arity == 1:
                    aggregate += paramDict[key].dot(alphaDict[tail[0]])
                elif arity == 2:
                    x1_alpha = alphaDict[tail[0]]
                    x2_alpha = alphaDict[tail[1]]
                    result = np.tensordot(paramDict[key], x1_alpha, axes=([1], [0]))
                    result = np.tensordot(result, x2_alpha, axes=([1], [0]))
                    aggregate += result
                else:
                    sys.stderr.write("Arity > 2! Cannot compute alpha terms\n")
        alphaDict[node.id] = aggregate
    return alphaDict

def computeOutside(hg, paramDict, rank, alphaDict):
    PDKeys = [pp for pp in paramDict if pp != "Pi"] #filtering sentence start parameters
    betaDict = {}
    for key in alphaDict: #initialize the beta terms to 0
        betaDict[key] = np.zeros((rank))
    betaDict[hg.nodes_[-1].id] = paramDict["Pi"] #initialization
    for q_node in reversed(hg.nodes_): #iterate through nodes in reverse topological order
        x_beta = betaDict[q_node.id]
        for inEdgeID in q_node.in_edges_: #for each incoming hyperedge
            src_rule = hg.edges_[inEdgeID].rule
            tail = hg.edges_[inEdgeID].tailNodes
            matching_rules = [rule for rule in PDKeys if rule.split(' ||| ')[1] == src_rule]
            for key in matching_rules: #for each rule where src RHS == hyper edge src RHS
                arity = checkArity(key)
                assert len(tail) == arity #skip if arity == 0 (pre-term)
                if arity == 1: #then just add the result to tailnode
                    betaDict[tail[0]] += x_beta.dot(paramDict[key])
                elif arity == 2:
                    result = np.tensordot(paramDict[key], x_beta, axes=([0], [0]))                
                    x_alpha_right = alphaDict[tail[1]]
                    betaDict[tail[0]] += np.tensordot(result, x_alpha_right, axes=([1], [0]))
                    x_alpha_left = alphaDict[tail[0]]
                    betaDict[tail[1]] += np.tensordot(result, x_alpha_left, axes=([0], [0]))                
    return betaDict
                
def insideOutside(hg, paramDict):
    rank = len(paramDict[paramDict.keys()[0]])
    alpha = computeInside(hg, paramDict, rank)
    beta = computeOutside(hg, paramDict, rank, alpha)
    marginals = {}
    for nodeID in alpha:
        marginals[nodeID] = np.dot(alpha[nodeID], beta[nodeID])
    return marginals
    

