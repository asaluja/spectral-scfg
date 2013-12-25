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
    return len(re.findall(r'\[([^]]+)\]', rule))

def computeInside(hg, paramDict, rank):
    alphaDict = {}
    for node in hg.nodes_: #nodes are added in topological order, so hg.nodes_ is guaranteed to be sorted
        aggregate = np.zeros((rank))
        for inEdgeID in node.in_edges_: #for each incoming hyperedge (i.e., rule that fits the span)
            src_rule = hg.edges_[inEdgeID].rule #src RHS string
            tail = hg.edges_[inEdgeID].tailNodes
            arity = checkArity(src_rule)
            assert len(tail) == arity
            key = ' ||| '.join([node.cat, src_rule])
            srcDict = paramDict[key] #dict of params with same src rule, each value is a tensor/matrix/vector
            for target_rule in srcDict: #key is target rule
                 if arity == 0: #pre-terminal
                     aggregate += srcDict[target_rule]
                 elif arity == 1:
                     aggregate += srcDict[target_rule].dot(alphaDict[tail[0]])
                 elif arity == 2:
                     x1_alpha = alphaDict[tail[0]]
                     x2_alpha = alphaDict[tail[1]]
                     result = np.tensordot(srcDict[target_rule], x1_alpha, axes=([1], [0]))
                     result = np.tensordot(result, x2_alpha, axes=([1], [0]))
                     aggregate += result
                 else:
                    sys.stderr.write("Arity > 2! Cannot compute alpha terms\n")
        alphaDict[node.id] = aggregate
    return alphaDict

def computeOutside(hg, paramDict, rank, alphaDict):
    betaDict = {}
    for key in alphaDict: #initialize the beta terms to 0
        betaDict[key] = np.zeros((rank))
    betaDict[hg.nodes_[-1].id] = paramDict["Pi"] #initialization
    for q_node in reversed(hg.nodes_): #iterate through nodes in reverse topological order
        x_beta = betaDict[q_node.id]
        for inEdgeID in q_node.in_edges_: #for each incoming hyperedge
            src_rule = hg.edges_[inEdgeID].rule #src RHS string
            tail = hg.edges_[inEdgeID].tailNodes
            arity = checkArity(src_rule)
            assert len(tail) == arity
            key = ' ||| '.join([q_node.cat, src_rule])
            srcDict = paramDict[key]
            for target_rule in srcDict: #for each rule where src RHS == hyper edge src RHS
                if arity == 1: #then just add the result to tailnode
                    betaDict[tail[0]] += x_beta.dot(srcDict[target_rule])
                elif arity == 2:
                    result = np.tensordot(srcDict[target_rule], x_beta, axes=([0], [0]))                
                    x_alpha_right = alphaDict[tail[1]]
                    betaDict[tail[0]] += np.tensordot(result, x_alpha_right, axes=([1], [0]))
                    x_alpha_left = alphaDict[tail[0]]
                    betaDict[tail[1]] += np.tensordot(result, x_alpha_left, axes=([0], [0]))                
    return betaDict
                
def insideOutside(hg, paramDict, rank):
    alpha = computeInside(hg, paramDict, rank)    #paramDict keys are 'LHS ||| src_RHS' format
    beta = computeOutside(hg, paramDict, rank, alpha)
    marginals = {}
    for nodeID in alpha:
        marginals[nodeID] = np.dot(alpha[nodeID], beta[nodeID])
        if marginals[nodeID] > 1 or marginals[nodeID] < 0:
            sys.stderr.write("Error! Marginal of span [%d,%d] outside of range: %.5g\n"%(hg.nodes_[nodeID].i, hg.nodes_[nodeID].j, marginals[nodeID]))
    return marginals
    

