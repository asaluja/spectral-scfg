#!/usr/bin/python -tt

import sys, commands, string, re

class trieNode:
    '''
    TrieNode class is the recursive data structure that trie
    relies on.  Each node contains a list of the rules that
    terminate at that node, and a dictionary that maps
    words (edges) to other trieNodes. 
    '''
    def __init__(self): #constructor for trieNode
        self.edges_ = {}
        self.rb_ = []
        
    def addEdge(self, item):
        if item not in self.edges_:
            self.edges_[item] = trieNode()

    def addRule(self, rule):
        self.rb_.append(rule)
        self.rb_ = list(set(self.rb_))

    def getRules(self):
        return self.rb_
            
    #returns a trieNode
    def Extend(self, word):
        if word in self.edges_:
            return self.edges_[word]
        else:
            return None

class trie:
    '''
    Trie class is a data structure used for decoding. 
    The trie class takes as input a list of all the rules in 
    the grammar, and constructs a trie where the edges correspond
    to words or NTs and each node contains a bin of all the rules
    that terminate at that node.  Given this data structure, we can 
    quickly retrieve the rules in the grammar that match a given 
    sub-span of a sentence.  
    To Do: 
    '''
    def __init__(self, rules): #constructor for trie
        self.root = trieNode()
        for rule in rules: #rules is a list of rules
            elements = rule.split(' ||| ')
            LHS = elements[0]
            RHS_src = elements[1]
            rule_src = self.formatRule(RHS_src.split()) #pass in RHS src
            curNode = self.root
            for item in rule_src: #add the rule to the trie
                curNode.addEdge(item)
                nextNode = curNode.Extend(item)
                curNode = nextNode
            curNode.addRule((LHS,RHS_src)) #once done, add rule to rule bin of current node   
    
    def formatRule(self, rule): #rule is the source phrase broken down into a list
        exp = re.compile(r'\[([^]]*)\]')
        rule_f = []
        for item in rule:
            if exp.match(item):
                rule_f.append('[X]')
            else:
                rule_f.append(item)
        return rule_f

    def hasRuleForSpan(self, start, end):
        return True

    def getRoot(self):
        return self.root
    
    '''for debugging'''
    def traverseTrie(self, pos, curNode=None):
        if pos == 0:
            curNode = self.root
        print "position %d rule bin: "%pos
        print curNode.getRules()
        for nextNode in curNode.edges_.keys():
            print "Traversing edge corresponding to %s"%(nextNode)
            self.traverseTrie(pos+1, curNode.Extend(nextNode))
            print "Returning to previous node"

class ActiveItem:
    '''
    ActiveItem is a data structure modeled on the ActiveItem struct in
    cdec's bottom_up_parser.cc.  
    '''
    def __init__(self, node, tailNodes=[]):
        self.gptr = node
        self.tailNodeVec = tailNodes

    def extendNonTerminal(self, node_idx):        
        ni = self.gptr.Extend("[X]")
        if ni is not None:
            tailNodes = self.tailNodeVec[:]
            tailNodes.append(node_idx)
            return ActiveItem(ni, tailNodes)
        else:
            return None
    
    def extendTerminal(self, word):
        ni = self.gptr.Extend(word)        
        if ni is not None:
            tail = self.tailNodeVec[:]
            return ActiveItem(ni, tail)
        else:
            return None

class HGEdge:
    '''
    This class represents the data structure for the edges in the generated
    hypergraph.  Each edge corresponds to the application of a rule.  
    Tail nodes is a vector of ints, where each int is the index of the nodes
    in the hypegraph nodes_ vector. The i and j variables correspond to the span
    of the head node, and are updated when we apply a rule from outside.  
    To Do: to make it more OOP-like, we might want to add a function that specifically
    modifies i and j.  
    '''
    def __init__(self, rule, tail, edgeID):
        self.rule = rule
        self.tailNodes = tail[:] #list of node IDs of the various tail nodes
        self.id = edgeID
        self.headNode = -1 #ID of head node

class HGNode:
    '''
    This class acts more as a struct, to encapsulate the various values
    (node ID, node category, ingoing and outgoing edges) associated
    with a given node in our hypergraph.  
    To Do: to make it more OOP-;like, we might want to add a function that
    specifically adds to outgoing or ingoing edges of a particular node. 
    '''
    def __init__(self, id_val, cat, start, end):
        self.id = id_val
        self.cat = cat
        self.in_edges_ = []
        self.out_edges_ = []
        self.i = start
        self.j = end

class HyperGraph:
    '''
    This class encapsulates the edges and nodes of the hypergraph 
    representation into two lists.  Remember that nodes correspond
    to NTs and have associated sub-spans, and edges correspond to
    the application of rules.  We also have utility functions to 
    add edges and nodes to the hypergraph, get a particular node 
    given an index (stored within the passive chart cells), connect
    an edge to its headNode, and a function to prune unreachable nodes. 
    '''
    def __init__(self):
        self.edges_ = []
        self.nodes_ = []

    def addEdge(self, rule, tailNodes):
        edge = HGEdge(rule, tailNodes, len(self.edges_))
        for nodeID in tailNodes: #append edge ID to outgoing edges of tail nodes
            self.nodes_[nodeID].out_edges_.append(edge.id)
        self.edges_.append(edge)
        return edge

    def addNode(self, cat, start, end):
        node = HGNode(len(self.nodes_), cat, start, end)
        self.nodes_.append(node)
        return node

    def getIthNode(self, idx):
        return self.nodes_[idx]

    def connectEdgeToHeadNode(self, edge, node):
        edge.headNode = node.id
        node.in_edges_.append(edge.id)
