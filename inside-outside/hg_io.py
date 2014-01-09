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
Update (Jan 2, 2013): added two different marginal computations: 
node marginal and edge marginal. 
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
                elif arity == 2: #this part of the code seems to take some time
                    x1_alpha = alphaDict[tail[0]]
                    x2_alpha = alphaDict[tail[1]]
                    result = np.tensordot(srcDict[target_rule], x2_alpha, axes=([2], [0]))
                    result = np.tensordot(result, x1_alpha, axes=([1], [0]))
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
                    result = np.tensordot(srcDict[target_rule], x_beta, axes=([0], [0])) #because we are defining axis, order of param and x_beta doesn't matter
                    x_alpha_right = alphaDict[tail[1]]
                    betaDict[tail[0]] += np.tensordot(result, x_alpha_right, axes=([1], [0]))
                    x_alpha_left = alphaDict[tail[0]]
                    betaDict[tail[1]] += np.tensordot(result, x_alpha_left, axes=([0], [0]))                
    return betaDict

def decorateSrcRule(hg, inEdgeID):
    expr = re.compile(r'\[([^]]*)\]')
    rule = hg.edges_[inEdgeID].rule
    tail = hg.edges_[inEdgeID].tailNodes[:]
    rule_decorated = []
    for item in rule.split():
        if expr.match(item): #NT, we need to decorate with its span
            child = hg.nodes_[tail.pop(0)]
            NT = child.cat[:-1] + "_%d_%d]"%(child.i,child.j)
            rule_decorated.append(NT)
        else:
            rule_decorated.append(item)
    return ' '.join(rule_decorated)

def decorateTgtRule(rule):
    expr = re.compile(r'\[([^]]*)\]')
    rule_decorated = []
    for item in rule.split():
        if expr.match(item): #NT, we need to strip down
            NTIdx = int(item.split(',')[1][0])
            rule_decorated.append("[%d]"%NTIdx)
        else:
            rule_decorated.append(item)
    return ' '.join(rule_decorated)

'''
in nodeMarginals, the marginal is associated with an NT span. 
This corresponds exactly to a node. 
'''
def nodeMarginals(alpha, beta, hg, flipSign, words):
    marginals = {}
    for node in hg.nodes_:
        LHS = node.cat[:-1] + "_%d_%d]"%(node.i, node.j)
        marginal = np.dot(alpha[node.id], beta[node.id])        
        if marginal < 0:
            if flipSign:
                marginal = -marginal
            else:
                sys.stderr.write("Error! Marginal of span [%d,%d] outside of range: %.5g\n"%(hg.nodes_[nodeID].i, hg.nodes_[nodeID].j, marginals[nodeID]))
                return marginals
        for inEdgeID in node.in_edges_: #the same marginal is used for all incoming edges (rules)
            key = ' ||| '.join([node.cat, hg.edges_[inEdgeID].rule])
            src_decorated = decorateSrcRule(hg, inEdgeID)
            for target_rule in paramDict[key]:
                src_tgt_decorated = "%s ||| %s"%(words[node.i], words[node.i]) if target_rule == "<unk>" else "%s ||| %s"%(src_decorated, decorateTgtRule(target_rule))
                lhs_src_tgt = ' ||| '.join([LHS, src_tgt_decorated])
                marginals[lhs_src_tgt] = marginal #associate the marginal to each rule
    return marginals

'''
in edge marginals, the marginal is associated with a particular
rule or edge.  We loop through all edges, obtain the beta vector
of the head node and the alpha vectors of the tail nodes, and
multiply them through the parameter associated with the rule to
get the edge marginal.  
'''
def edgeMarginals(alpha, beta, hg, flipSign, paramDict, words):
    marginals = {}
    for edge in hg.edges_:
        head = hg.nodes_[edge.headNode]
        tail = edge.tailNodes        
        beta_head = beta[edge.headNode]
        LHS = head.cat[:-1] + "_%d_%d]"%(head.i, head.j)
        key = ' ||| '.join([head.cat, edge.rule])
        src_decorated = decorateSrcRule(hg, edge.id)
        arity = checkArity(edge.rule)
        assert arity == len(tail)
        for target_rule in paramDict[key]:
            #src_tgt_decorated = "%s ||| %s"%(words[head.i], words[head.i]) if target_rule == "<unk>" else "%s ||| %s"%(src_decorated, decorateTgtRule(target_rule))
            src_tgt_decorated = "<unk> ||| %s"%(words[head.i]) if target_rule == "<unk>" else "%s ||| %s"%(src_decorated, decorateTgtRule(target_rule))
            lhs_src_tgt = ' ||| '.join([LHS, src_tgt_decorated])
            marginal = 0
            if arity == 0:
                marginal = beta_head.dot(paramDict[key][target_rule])
            elif arity == 1:
                marginal = beta_head.dot(paramDict[key][target_rule]).dot(alpha[tail[0]])
            elif arity == 2:
                result = np.tensordot(paramDict[key][target_rule], beta_head, axes=([0], [0])) 
                marginal = alpha[tail[0]].dot(result).dot(alpha[tail[1]])
            else:
                sys.stderr.write("Arity > 2! Cannot compute marginals for rule %s\n"%lhs_src_tgt)
            if marginal < 0:
                if flipSign:
                    marginal = -marginal
                else:
                    sys.stderr.write("Error! Marginal of span [%d,%d] outside of range: %.5g\n"%(head.i, head.j, marginal))
            if marginal > 0:
                marginals[lhs_src_tgt] = marginal
    return marginals
                
def insideOutside(hg, paramDict, rank, words, flipSign, nodeMarginal):
    alpha = computeInside(hg, paramDict, rank)    #paramDict keys are 'LHS ||| src_RHS' format
    beta = computeOutside(hg, paramDict, rank, alpha)
    marginals = nodeMarginals(alpha, beta, hg, flipSign, words) if nodeMarginal else edgeMarginals(alpha, beta, hg, flipSign, paramDict, words)
    return marginals #key of marginals is a decorated entire rule LHS ||| src RHS ||| tgt RHS
    

