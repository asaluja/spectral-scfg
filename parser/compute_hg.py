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
Update: incorporated python's multiprocessing setup, much more efficient than the
previous multiprocessing setup.  
Usage: python intersect_scfg.py (-d/-f/-n) SpectralParams Rank InputFile NumProcesses outDir/
Update (May 20,m 2014): renamed the directory from 'inside-outside' to 'parser', and also 
renamed this file from 'intersect_scfg.py' to 'compute_hg.py'.  Changed writing out of the marginals,
instead of log marginals previously we have raw marginals.  Also, there is now an option to write
out multiple spectral marginals (normalized by source or target). 
Lastly, included an additional option flag to generate heat maps of the parse chart. 
Usage: python compute_hg.py (-d/-f/-n/-m/-s/-t) params rank input_sentences numProc outDir/
'''

import sys, commands, string, time, gzip, cPickle, re, getopt, math, hg_io
import multiprocessing as mp
import pylab as plt
import matplotlib as mpl
import numpy as np
from trie import trie, ActiveItem, HyperGraph

(opts, args) = getopt.getopt(sys.argv[1:], 'dnfmstx')
optsDict = {}
for opt in opts:
    if opt[0] == '-d': #print info about each node in the hypergraph 
        optsDict["debug"] = 1
    elif opt[0] == '-n': #by default, we print edge marginal
        optsDict["nodeMarginal"] = 1 #if true, we print out heat maps
    elif opt[0] == '-f': #if marginal is < 0, we flip the sign
        optsDict["flipSign"] = 1
    elif opt[0] == '-m': #MLE 
        optsDict["MLE"] = 1
    elif opt[0] == '-s': #source norm
        optsDict["sourceNorm"] = 1
    elif opt[0] == '-t': #target norm
        optsDict["targetNorm"] = 1
    elif opt[0] == '-x': #only write out marginals for non-lexical rules in source
        optsDict["onlyXX"] = 1

params_fh = open(args[0], 'rb')
paramDict = cPickle.load(params_fh) #key is 'LHS ||| src RHS'
grammar_rules = [rule for rule in paramDict.keys() if rule != "Pi"] #Pi contains the start of sentence params
grammarTrie = trie(grammar_rules) 
rank = int(args[1])
inputFile = open(args[2], 'r').readlines()
numProcesses = int(args[3])
outDir = args[4]

'''
declaration of list that maintains which sentences have failed across all processes
'''
def init(fs, fls):
    global failed_sentences, flipped_sentences
    failed_sentences = fs
    flipped_sentences = fls

def main():
    failed_sentences = mp.Manager().list()
    flipped_sentences = mp.Manager().list()
    pool = mp.Pool(processes=numProcesses, initializer=init, initargs=(failed_sentences, flipped_sentences))
    for i, line in enumerate(inputFile):
        outFile = "%s/chart.%d.pdf"%(outDir, i) if "nodeMarginal" in optsDict else "%s/grammar.%d.gz"%(outDir, i)
        #parse(line.strip().split(), outFile, i)
        pool.apply_async(parse, (line.strip().split(), outFile, i))
    pool.close()
    pool.join()
    print "number of failed sentences: %d"%(len(failed_sentences))
    print "number of flipped sentences: %d"%(len(flipped_sentences))

'''
main function for bottom-up parser with Earley-style rules. 
The active chart is first seeded with pointers to the root
node of a source rules trie. Then, in a bottom-up manner, 
we advance the dots for each cell item, and then convert completed
rules in a cell to the passive chart, or deal with NTs in active
items just proved.  At the end, we look at the passive items in
the cell corresponding to the sentence to see if [S] is there. 
'''
def parse(words, outFile, lineNum):
    start = time.clock()
    N = len(words)
    goal_idx = False
    hg = HyperGraph()
    active = {}
    passive = {}
    nodemap = {}
    seedActiveChart(N, active)
    for l in range(1, N+1): #length
        for i in range(0, N+1-l): #left index of span            
            j = i + l #right index of span
            advanceDotsForAllItemsInCell(i, j, words, active, passive)
            cell = active[(i,j)][:] if (i,j) in active else [] #list of active items
            for activeItem in cell:
                rules = activeItem.srcTrie.getRules()
                for rule in rules:
                    applyRule(i, j, rule, activeItem.tailNodeVec, hg, nodemap, passive)
            if j < N: #the below function includes NTs that were just proved into new binaries, which is unnecessary for the end token
                extendActiveItems(i, i, j, active, passive) #dealing with NTs that were just proved
        if (0,N) in passive: #we have spanned the entire input sentence
            passiveItems = passive[(0,N)] #list of indices
            if len(passiveItems) > 0: #we have at least one node that covers the entire sentence
                goal_idx = True                
    parseTime = time.clock() - start
    if goal_idx: #i.e., if we have created at least 1 node in the HG corresponding to goal        
        print "SUCCESS; length: %d words, time taken: %.2f sec, sentence: %s"%(len(words), parseTime, ' '.join(words))
    else:
        print "FAIL; length: %d words, time taken: %.2f sec"%(len(words), parseTime)
        failed_sentences.append(lineNum)
    if "debug" in optsDict:
        printHGDebug(hg, lineNum)
    elif goal_idx:
        computeMarginals(hg, words, outFile, lineNum)
    sys.stdout.flush()

def computeMarginals(hg, words, outFile, lineNum):
    flipSign = "flipSign" in optsDict
    nodeMarginal = "nodeMarginal" in optsDict
    sourceNorm = "sourceNorm" in optsDict
    targetNorm = "targetNorm" in optsDict
    MLE = "MLE" in optsDict
    XX = "onlyXX" in optsDict
    rankToUse = 0 if MLE else rank
    start = time.clock()
    marginals, flipped = hg_io.insideOutside(hg, paramDict, rankToUse, words, flipSign, nodeMarginal)
    if marginals is not None:
        if flipped:
            flipped_sentences.append(lineNum)
        ioTime = time.clock() - start
        print "marginals computed over hypergraph. time taken: %.2f sec"%(ioTime)
        if nodeMarginal:
            printHeatMap(marginals, words, outFile)
        else:
            tgt_norm_marginals = None
            src_norm_marginals = None
            if targetNorm:
                tgt_norm_marginals = normalizeMarginalsByTarget(marginals)
            if sourceNorm:
                src_norm_marginals = normalizeMarginalsBySource(marginals)        
            fh = gzip.open(outFile, 'w')
            for key in marginals:
                featureStr = "%s ||| spectral=%.5g"%(key, marginals[key])
                if sourceNorm:
                    featureStr += " spectralEgivenF=%.5g"%src_norm_marginals[key]
                if targetNorm:
                    featureStr += " spectralFgivenE=%.5g"%tgt_norm_marginals[key]
                if XX: 
                    if checkNonLex(key): #if non-lexical, write out the features
                        fh.write("%s\n"%featureStr)
                    else: #otherwise, just write out the rule
                        fh.write("%s\n"%(key))
                else:
                    fh.write("%s\n"%featureStr)
            fh.write("[S] ||| [X_0_%d] ||| [1] ||| 0\n"%len(words)) #top level rule
            fh.close()

def checkNonLex(rule):
    elements = rule.split(' ||| ')
    numNTs = len(re.findall(r'\[([^]]+)\]', elements[1]))
    numWords = len(elements[1].split())
    if numNTs == numWords:
        return True
    else:
        return False

def printHeatMap(marginals, words, outFile):
    N = len(words)
    words_uni = [i.decode('UTF-8') for i in words]
    heatmap = np.zeros((N+1, N+1))
    for chart in marginals:
        heatmap[chart[0], chart[1]] = math.log(marginals[chart])
    fig, ax = plt.subplots()    
    mask = np.tri(heatmap.shape[0], k=0)
    heatmap = np.ma.array(heatmap, mask=mask)
    cmap = plt.cm.get_cmap('RdBu')
    cmap.set_bad('w')
    im = ax.pcolor(heatmap, cmap=cmap, alpha=0.8)
    font = mpl.font_manager.FontProperties(fname='/usr0/home/avneesh/spectral-scfg/data/wqy-microhei.ttf')
    ax.grid(True)
    ax.set_ylim([0,N])
    ax.invert_yaxis()
    ax.set_yticks(np.arange(heatmap.shape[1]-1)+0.5, minor=False)
    ax.set_yticklabels(words_uni, minor=False, fontproperties=font)
    ax.set_xticks(np.arange(heatmap.shape[0])+0.5, minor=True)
    ax.set_xticklabels(np.arange(heatmap.shape[0]), minor=True)
    ax.set_xticks([])
    cbar = fig.colorbar(im, use_gridspec=True)
    cbar.set_label('ln(sum)')
    ax.set_xlabel('Span End')
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    plt.ylabel('Span starting at word: ')
    plt.tight_layout()
    #ax.set_title('CKY Heat Map: Node Marginals')
    fig.savefig(outFile)    


def normalizeMarginalsByTarget(margDict):
    sortMargDict = {}
    for key in margDict:
        elements = key.split(' ||| ')
        target = ' ||| '.join([elements[0], elements[2]]) #target is actually LHS + RHS target
        src_marginal = (elements[1], margDict[key])
        if target in sortMargDict:
            sortMargDict[target].append(src_marginal)
        else:
            sortMargDict[target] = [src_marginal]
    normMargDict = {}
    for key in sortMargDict:
        normalizer = sum([src_marginal[1] for src_marginal in sortMargDict[key]])
        LHS = key.split(' ||| ')[0]
        target = key.split(' ||| ')[1]
        for src_marginal in sortMargDict[key]: #put it back in normal order            
            format_key = "%s ||| %s ||| %s"%(LHS, src_marginal[0], target)
            normMargDict[format_key] = 0.999 if src_marginal[1]/normalizer == 1 else src_marginal[1]/normalizer
    return normMargDict 

def normalizeMarginalsBySource(margDict):
    sortMargDict = {}
    for key in margDict:
        elements = key.split(' ||| ')
        source = ' ||| '.join([elements[0], elements[1]])
        tgt_marginal = (elements[2], margDict[key])
        if source in sortMargDict:
            sortMargDict[source].append(tgt_marginal)
        else:
            sortMargDict[source] = [tgt_marginal]
    normMargDict = {}
    for key in sortMargDict:
        normalizer = sum([tgt_marginal[1] for tgt_marginal in sortMargDict[key]])
        for tgt_marginal in sortMargDict[key]:
            format_key = "%s ||| %s"%(key, tgt_marginal[0])
            normMargDict[format_key] = tgt_marginal[1]/normalizer
    return normMargDict

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
Function called before the sentence is parsed;
places a pointer to the source rules trie root
along the diagonal of the active chart. 
'''
def seedActiveChart(N, active):
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
def advanceDotsForAllItemsInCell(start, end, words, active, passive):
    for k in range(start+1, end):
        extendActiveItems(start, k, end, active, passive)        
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
def extendActiveItems(start, split, end, active, passive):
    icell = active[(start, split)] if (start,split) in active else []
    idxs = passive[(split, end)] if (split,end) in passive else []
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
def applyRule(start, end, rule, tailNodes, hg, nodemap, passive):
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
            

    
    
    
    
    
    
    
    
    
    
        
            

            
    
    
    
    
    
    

