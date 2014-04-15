#!/usr/bin/python -tt

'''
File: spectral_context.py
Date: January 9, 2014
Description: 
'''

import sys, commands, string, os, gzip, re, itertools, operator
import scipy.stats as stats
from collections import Counter

def extractPreTerminals(filehandle, N):
    rules = [[] for x in xrange(0, N)]
    expr = re.compile(r'\[([^]]*)\]')
    for rule in filehandle:
        elements = rule.strip().split(' ||| ')
        src = elements[1]
        if len(src.split()) == 1 and not expr.match(src): #one-word pre-term and not an NT
            pos = int(elements[0].split('_')[1])
            marginal = float(elements[3].split('=')[1])
            rules[pos].append((elements[2], marginal))
    return rules

def storeSpectralByContext(sentence, rules, contextDict):
    for i in xrange(0, len(rules)): #go through every rule position
        left_context = "<s>" if i == 0 else sentence[i-1]
        right_context = "</s>" if i == len(rules) - 1 else sentence[i+1]
        word = sentence[i]
        contexts = contextDict[word] if word in contextDict else {}
        sorted_marginals = sorted(rules[i], key=operator.itemgetter(1), reverse=True) #list of tuples, each entry of the form (tgtRule, marginal)
        key = "%s %s %s"%(left_context, word, right_context)
        contexts.setdefault(key, []).append(sorted_marginals) #adds sorted list 
        contextDict[word] = contexts

'''
each entry in rankings is a list of tuples (target word, marginal) that occur for this context
function is only called when len(rankings) > 1, i.e., the trigram has been seen more than once
'''
def computeTauForAllPairs(rankings): 
    pairDict = {}
    for pair in list(itertools.combinations(range(len(rankings)), 2)): #loop through all index pairs
        rank1 = [rank[0] for rank in rankings[pair[0]]] #extract just the sorted words
        rank2 = [rank[0] for rank in rankings[pair[1]]]
        assert len(rank1) == len(rank2) #the number of rules over different targets should be the same across contexts!
        if len(rank1) == 1: #only one target translation
            pairDict[pair] = 1 if rank1[0] == rank2[0] else -1 #1 if they're the same, -1 if they're different
        else:
            tau, pval = stats.kendalltau(rank1, rank2) #higher tau means strong agreement
            pairDict[pair] = tau
    #figure out which pairs have the highest and lowest tau    
    sorted_pairs = sorted(pairDict.iteritems(), key=operator.itemgetter(1)) #sort pairs in ascending order based on tau; higher tau means more agreement
    return sorted_pairs[0][1] #return actual tau value

'''
each entry in rankings is a tuple (trigram, most frequent ranking given trigram) for a given word
'''
def computeTauAcrossContexts(rankings, word):
    pairDict = {}
    for pair in list(itertools.combinations(range(len(rankings)), 2)): #loop through all pairs of contexts
        rank1 = rankings[pair[0]][1] #list over target side words
        rank2 = rankings[pair[1]][1]        
        key1 = rankings[pair[0]][0]
        key2 = rankings[pair[1]][0]
        assert len(rank1) == len(rank2)
        assert len(rank1) > 0 #needs to have at least 1 target side translation
        if len(rank1) == 1: #only one target translation
            pairDict[(key1, key2)] = 1 if rank1[0] == rank2[0] else -1 #1 if they're the same, -1 if theyr'e different
        else:
            tau, pval = stats.kendalltau(rank1, rank2) #higher tau means strong agreement
            pairDict[(key1, key2)] = tau            
    sorted_pairs = sorted(pairDict.iteritems(), key=operator.itemgetter(1)) #returns ascending list of tuples ((trigram 1, trigram 2), tau)
    return sorted_pairs[0] #returns most divergent context pair, and its tau

def analyzeSpectralContext(contextDict):
    withinContextVar = {}
    acrossContexts = {}
    for word in contextDict:
        withinContextVar[word] = {}
        across_context_rankings = []
        for context in contextDict[word]: #context is a trigram string
            assert len(contextDict[word][context]) > 0
            if len(contextDict[word][context]) > 1: #meaning we have encountered this trigram more than once in the data
                lowest_tau = computeTauForAllPairs(contextDict[word][context])
                if lowest_tau < 1: 
                    print "context '%s': observed different rankings over target phrases in the same context"%(context)
                sorted_word_rankings = [tuple([wordSpectralPair[0] for wordSpectralPair in ranking]) for ranking in contextDict[word][context]] 
                most_freq_ranking = list(Counter(sorted_word_rankings).most_common()[0][0]) #find the mode and use that as the representative ranking
                #print "most frequent ranking: "
                #print most_freq_ranking
                #check and see what most_freq_ranking is --> should be the same type as sorted_word_ranking
                across_context_rankings.append((context, most_freq_ranking))               
            else:                
                RuleMarginalPairList = contextDict[word][context][0]
                sorted_word_ranking = [RuleMarginalPair[0] for RuleMarginalPair in RuleMarginalPairList]
                #print "sorted word ranking: "
                #print sorted_word_ranking
                across_context_rankings.append((context, sorted_word_ranking))
        if len(across_context_rankings) > 1: #seen at least two different contexts for this word; otherwise the analysis is irrelevant
            acrossContexts[word] = computeTauAcrossContexts(across_context_rankings, word) #provides most divergent context pair, and the tau
    for word in acrossContexts:
        contextPair, div = acrossContexts[word]
        key1, key2 = contextPair
        print "word '%s' ||| context pair '%s' and '%s' ||| Kendall's Tau: %.3f"%(word, key1, key2, div)

def printContextDict(contextDict):
    for word in contextDict:
        for context in contextDict[word]:
            print "context: %s"%(context)
            print contextDict[word][context]

def main():
    mingrammar_loc = sys.argv[1]
    contextDict = {}
    for i,line in enumerate(sys.stdin): #loop through all input sentences
        sentence = line.strip().split()
        minrule_fh = gzip.open(mingrammar_loc + "grammar.%d.gz"%i, 'rb')
        rules_by_pos = extractPreTerminals(minrule_fh, len(sentence))
        storeSpectralByContext(sentence, rules_by_pos, contextDict)
    #printContextDict(contextDict)
    analyzeSpectralContext(contextDict)

if __name__ == "__main__":
    main()
