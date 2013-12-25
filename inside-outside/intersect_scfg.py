#!/usr/bin/python -tt
'''
File: intersect_scfg.py
Date: December 13, 2013
Description: modified version of bottom-up parser
where we use a variant of CKY+, leveraging a trie 
data structure of the grammar to come up with a hypergraph
representation of the composition of the input sentence
FST and the grammar. The hypergraph is then converted
to a grammar where NTs are indexed/decorated with their span.
The grammar is written out in .gz files separated by sentence. 
Downstream, io.py must read in and process the per-sentence grammar
into a chart, after which the alpha and beta terms can be computed easily. 
arg1: dictionary of parameters (output of feature extraction step)
stdin: tokenized sentences
'''

import sys, commands, string, cPickle, re, hg_io
from trie import trie, ActiveItem, HyperGraph

params_fh = open(sys.argv[1], 'rb')
paramDict = cPickle.load(params_fh) #key is 'LHS ||| src RHS'
rank = int(sys.argv[2])
passive = {}
active = {}
nodemap = {}

def main():
    grammar_rules = [rule for rule in paramDict.keys() if rule != "Pi"] #Pi contains the start of sentence params
    grammarTrie = trie(grammar_rules) 
    for line in sys.stdin:
        print "bottom-up parse for input sentence: "
        print line.strip()
        words = line.strip().split()
        goal_nodes = []
        hg = parse(words, grammarTrie, goal_nodes)
        active.clear()
        passive.clear()
        nodemap.clear()
        if len(goal_nodes) > 0: #i.e., if we have created at least 1 node in the HG corresponding to goal
            print "parsing success; length of sentence: %d"%(len(words))
            #print "(Nodes/Edges): %d / %d"%(len(hg.nodes_), len(hg.edges_))
            marginals = hg_io.insideOutside(hg, paramDict, rank)
            print "marginals computed over hypergraph"
            convertHGToRules(hg, marginals)
        else:
            print "parsing fail; length of sentence: %d"%(len(words))
        sys.stdout.flush()
 
def convertHGToRules(hg, marginals):
    expr = re.compile(r'\[([^]]*)\]')
    for node in hg.nodes_:        
        LHS = node.cat[:-1] + ",%d-%d]"%(node.i, node.j)
        marginal = marginals[node.id] #marginals defined over span
        for inEdgeID in node.in_edges_:
            rule = hg.edges_[inEdgeID].rule
            tail = hg.edges_[inEdgeID].tailNodes[:]
            rule_decorated = []
            for item in rule.split():
                if expr.match(item): #NT, we need to decorate with its span
                    child = hg.nodes_[tail.pop(0)]
                    NT = child.cat[:-1] + ",%d-%d]"%(child.i,child.j)
                    rule_decorated.append(NT)
                else:
                    rule_decorated.append(item)
            print "%s ||| %s ||| %.5g"%(LHS, ' '.join(rule_decorated), marginal)

def parse(words, grammarTrie, goal_idx):
    N = len(words)
    hg = HyperGraph()
    seedActiveChart(grammarTrie, N)
    for l in range(1, N+1): #length
        for i in range(0, N+1-l): #left index of span            
            j = i + l #right index of span
            advanceDotsForAllItemsInCell(i, j, words)
            cell = active[(i,j)][:] if (i,j) in active else [] #list of active items
            for activeItem in cell:
                rules = activeItem.gptr.getRules()
                #if we are at (0,N), then we should only apply rules that have 'S' as LHS; otherwise, we should only apply rules that have 'X' as LHS
                rules = [rule for rule in rules if rule[0] == '[S]'] if (i == 0) and (j == N) else [rule for rule in rules if rule[0] == '[X]'] 
                if len(rules) > 0:
                    applyRules(i, j, rules, activeItem.tailNodeVec, hg)
            if j < N: #the below function includes NTs that were just proved into new binaries, which is unnecessary for the end token
                extendActiveItems(i, i, j) #dealing with NTs that were just proved
        if (0,N) in passive: #we have spanned the entire input sentence
            passiveItems = passive[(0,N)] #list of indices
            for idx in passiveItems:
                node = hg.getIthNode(idx)
                if node.cat == '[S]': #then goal state has been found
                    goal_idx.append(node.id)
    return hg

'''
places a pointer to the source rules trie root along the diagonal of the active chart
'''
def seedActiveChart(grammarTrie, N):
    global active
    for i in range(0, N): #note: for now, we don't test hasRuleForSpan        
        active.setdefault((i,i), []).append(ActiveItem(grammarTrie.getRoot())) #add the root of the trie

def advanceDotsForAllItemsInCell(start, end, words):
    for k in range(start+1, end):
        extendActiveItems(start, k, end)        
    ec = active[(start,end-1)] if (start,end-1) in active else []
    word = words[end-1]
    for actItem in ec:
        ai = actItem.extendTerminal(word)
        if ai is not None:
            active.setdefault((start,end), []).append(ai)

def extendActiveItems(start, split, end):
    global active
    icell, idxs = [], []
    if (start,split) in active:
        icell = active[(start, split)]
    if (split,end) in passive:
        idxs = passive[(split, end)]
    for actItem in icell:
        for idx in idxs:
            ai = actItem.extendNonTerminal(idx) 
            if ai is not None:
                active.setdefault((start,end), []).append(ai)

def applyRules(start, end, rules, tailNodes, hg):
    for rule in rules:
        applyRule(start, end, rule, tailNodes, hg)

def applyRule(start, end, rule, tailNodes, hg):
    global nodemap, passive
    edge = hg.addEdge(rule[1], tailNodes) #rule[1] is src RHS of rule
    node = None
    cat2NodeMap = {}
    if (start,end) in nodemap:
        cat2NodeMap = nodemap[(start,end)]    
    LHS = rule[0]
    if LHS in cat2NodeMap: #LHS is either [X] or [S] --> test if this ever fires?
        node = hg.nodes_[cat2NodeMap[LHS]]
    else:
        node = hg.addNode(LHS, start, end)
        cat2NodeMap[LHS] = node.id
        nodemap[(start,end)] = cat2NodeMap
        passive.setdefault((start,end), []).append(node.id)
    hg.connectEdgeToHeadNode(edge, node)

if __name__ == "__main__":
    main()
            

    
    
    
    
    
    
    
    
    
    
        
            

            
    
    
    
    
    
    

