#!/usr/bin/python -tt

import sys, commands, string, cPickle
from trie import trie

param_fh = open(sys.argv[1], 'rb')
paramDict = cPickle.load(param_fh)
srcPhrases = [key[0] for key in paramDict.keys() if key[0] != 'S']
grammarTrie = trie(srcPhrases)
grammarTrie.traverseTrie(0)
