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
Usage: python intersect_scfg.py (-d/f/n) SpectralParams SentencesToDecode NumPartitions Partition Rank OutDir
'''

import sys, commands, string, time, gzip, cPickle, re, getopt, math, hg_io
from trie import trie, ActiveItem, HyperGraph

(opts, args) = getopt.getopt(sys.argv[1:], 'dnf')
optsDict = {}
for opt in opts:
    if opt[0] == '-d': #print info about each node in the hypergraph 
        optsDict["debug"] = 1
    elif opt[0] == '-n': #by default, we print edge marginal
        optsDict["nodeMarginal"] = 1
    elif opt[0] == '-f': #if marginal is < 0, we flip the sign
        optsDict["flipSign"] = 1

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
    failed_sentences = []
    for lineNum in range(start,end):
        line = inputFile[lineNum]
        words = line.strip().split()
        start = time.clock()
        hg, goalReached = parse(words, grammarTrie)
        parseTime = time.clock() - start
        active.clear()
        passive.clear()
        nodemap.clear()
        if goalReached: #i.e., if we have created at least 1 node in the HG corresponding to goal
            print "SUCCESS; length: %d words, time taken: %.2f sec, sentence: %s"%(len(words), parseTime, line.strip())
            if "debug" in optsDict:
                printHGDebug(hg, lineNum)
            else:
                flipSign = "flipSign" in optsDict
                nodeMarginal = "nodeMarginal" in optsDict
                start = time.clock()
                marginals = hg_io.insideOutside(hg, paramDict, rank, words, flipSign, nodeMarginal)
                ioTime = time.clock() - start
                print "marginals computed over hypergraph. time taken: %.2f sec"%(ioTime)
                fh = gzip.open("%s/grammar.%d.gz"%(outDir, lineNum), 'w')
                for key in marginals:
                    fh.write("%s ||| spectral=%.3f\n"%(key, math.log(marginals[key])))
                fh.write("[S] ||| [X_0_%d] ||| [1] ||| 0\n"%len(words)) #top level rule
                fh.close()
        else:
            print "FAIL; length: %d words, time taken: %.2f sec"%(len(words), parseTime)
            failed_sentences.append(lineNum)
            if "debug" in optsDict:
                printHGDebug(hg, lineNum)
        sys.stdout.flush()
    print "number of failed sentences: %d"%(len(failed_sentences))

'''
hypergraph print function for debugging purposes
'''
def printHGDebug(hg, line_num):
    fh = open("%s/%d"%(outDir, line_num), 'w')
    print >> fh, "(Nodes/Edges): %d / %d"%(len(hg.nodes_), len(hg.edges_))
    for node in hg.nodes_: 
        print >> fh, "Node ID: %d"%(node.id)
        LHS = node.cat[:-1] + "_%d_%d]"%(node.i, node.j)
        for inEdgeID in node.in_edges_:
            srcRuleDecorated = hg_io.decorateSrcRule(hg, inEdgeID)
            print >> fh, "%s ||| %s"%(LHS, srcRuleDecorated)
    fh.close()

'''
main function for bottom-up parser with Earley-style rules. 
The active chart is first seeded with pointers to the root
node of a source rules trie. Then, in a bottom-up manner, 
we advance the dots for each cell item, and then convert completed
rules in a cell to the passive chart, or deal with NTs in active
items just proved.  At the end, we look at the passive items in
the cell corresponding to the sentence to see if [S] is there. 
'''
def parse(words, grammarTrie):
    N = len(words)
    goal_idx = False
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
                #rules = [rule for rule in rules if rule[0] == '[S]'] if (i == 0) and (j == N) else [rule for rule in rules if rule[0] == '[X]'] 
                for rule in rules:
                    applyRule(i, j, rule, activeItem.tailNodeVec, hg)
            if j < N: #the below function includes NTs that were just proved into new binaries, which is unnecessary for the end token
                extendActiveItems(i, i, j) #dealing with NTs that were just proved
        if (0,N) in passive: #we have spanned the entire input sentence
            passiveItems = passive[(0,N)] #list of indices
            if len(passiveItems) > 0: #we have at least one node that covers the entire sentence
                goal_idx = True
    return (hg, goal_idx)

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
        if end-start == 1: #OOV handling
            if ai is None:
                active.setdefault((start,end), []).append(actItem.extendOOV())                
            else: #check if active item has any rules in its bin
                if len(ai.srcTrie.getRules()) == 0: #handles the case where rule starts with OOV word, but no rule that actually covers OOV word
                    active.setdefault((start,end), []).append(actItem.extendOOV())
'''
function that extends active items over non-terminals. 
'''
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

'''
Given a rule, does the necessary book-keeping to 
convert that rule to the passive chart, and adds the 
appropriate nodes and edges to the hypergraph. 
'''
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
            

    
    
    
    
    
    
    
    
    
    
        
            

            
    
    
    
    
    
    

