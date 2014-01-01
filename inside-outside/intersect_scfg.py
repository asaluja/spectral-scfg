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
Update: incorporated simple multicore setup.  Basically, we divide the input
test corpus into partitions or chunks, and each process handles one partition. 
Usage: python intersect_scfg.py (-d) SpectralParams SentencesToDecode NumPartitions Partition Rank OutDir
'''

import sys, commands, string, cPickle, re, getopt, math, hg_io
from trie import trie, ActiveItem, HyperGraph

(opts, args) = getopt.getopt(sys.argv[1:], 'd')
debug=False
for opt in opts:
    if opt[0] == '-d':
        debug=True

params_fh = open(args[0], 'rb')
paramDict = cPickle.load(params_fh) #key is 'LHS ||| src RHS'
inputFile = open(args[1], 'r').readlines()
numPartitions = int(args[2])
partition = int(args[3])
rank = int(args[4])
outDir = args[5]
passive = {}
active = {}
nodemap = {}

def main():
    grammar_rules = [rule for rule in paramDict.keys() if rule != "Pi"] #Pi contains the start of sentence params
    grammarTrie = trie(grammar_rules) 
    numSentences = len(inputFile)
    sentencesPerChunk = numSentences/numPartitions #int/int means sentencesPerChunk will be truncated
    start = partition*sentencesPerChunk
    end = (partition+1)*sentencesPerChunk if partition < numPartitions - 1 else numSentences
    for lineNum in range(start,end):
        line = inputFile[lineNum]
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
            if debug:
                printHGDebug(hg, lineNum)
            else:
                marginals = hg_io.insideOutside(hg, paramDict, rank)
                print "marginals computed over hypergraph"
                convertHGToRules(hg, marginals, lineNum, words)
        else:
            print "parsing fail; length of sentence: %d"%(len(words))
        sys.stdout.flush()

def printHGDebug(hg, line_num):
    print "(Nodes/Edges): %d / %d"%(len(hg.nodes_), len(hg.edges_))
    fh = open("%s/%d"%(outDir, line_num), 'w')
    for node in hg.nodes_: 
        print >> fh, "Node ID: %d"%(node.id)
        LHS = node.cat[:-1] + ",%d-%d]"%(node.i, node.j)
        for inEdgeID in node.in_edges_:
            rule_decorated = decorateRule(hg, inEdgeID)
            print >> fh, "%s ||| %s"%(LHS, rule_decorated)
    fh.close()

def convertHGToRules(hg, marginals, line_num, words):
    fh = open("%s/%d"%(outDir, line_num), 'w')
    for node in hg.nodes_: 
        LHS = node.cat[:-1] + ",%d-%d]"%(node.i, node.j)
        marginal = marginals[node.id] #marginals defined over span
        if marginal > 0:            
            for inEdgeID in node.in_edges_:
                key = ' ||| '.join([node.cat, hg.edges_[inEdgeID].rule])
                src_decorated = decorateRule(hg, inEdgeID)
                for target_rule in paramDict[key]: 
                    src_tgt_decorated = "%s ||| %s"%(words[node.i], words[node.i]) if target_rule == "<unk>" else "%s ||| %s"%(src_decorated, target_rule)
                    print >> fh, "%s ||| %s ||| %.5g"%(LHS, src_tgt_decorated, math.log(marginal))
    fh.close()

def decorateRule(hg, inEdgeID):
    expr = re.compile(r'\[([^]]*)\]')
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
    return ' '.join(rule_decorated)

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
                rules = activeItem.srcTrie.getRules()
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
Function called before the sentence is parsed;
places a pointer to the source rules trie root
along the diagonal of the active chart. 
'''
def seedActiveChart(grammarTrie, N):
    global active
    for i in range(0, N): #note: for now, we don't test hasRuleForSpan        
        active.setdefault((i,i), []).append(ActiveItem(grammarTrie.getRoot())) #add the root of the trie

'''
Function that "advances the dot" (in a dotted rule)
on position to the right for all active items in the cell
defined by (start,end).  We first perform online binarization
by looping through all split points in the span and then see if
advancing the dot happened to cover a non-terminal (this is handled
in extendActiveItems).  We then check and see if advancing the dot 
happened to cover a new rule with the additional terminal.  
'''
def advanceDotsForAllItemsInCell(start, end, words):
    for k in range(start+1, end):
        extendActiveItems(start, k, end)        
    ec = active[(start,end-1)] if (start,end-1) in active else []
    word = words[end-1]
    for actItem in ec:
        ai = actItem.extendTerminal(word)
        if ai is not None:
            active.setdefault((start,end), []).append(ai)
        elif end-start == 1: #ai is None but we are covering an OOV word
            print "Extended OOV for span %d,%d"%(start,end)
            active.setdefault((start,end), []).append(actItem.extendOOV())

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
            

    
    
    
    
    
    
    
    
    
    
        
            

            
    
    
    
    
    
    

