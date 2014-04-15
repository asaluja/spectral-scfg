#!/usr/bin/python -tt

import sys, commands, string. gzip, getopt, cPickle, os

def updateCountsFromExample(stree):

def main():
    optsDict = {}
    (opts, args) = getopt.getopt(sys.argv[1:], 'n:')
    for opt in opts:
        if opt[0] == '-n': #number of iterations
            optsDict["iterations"] = int(opt[1])            
    minrule_grammars_loc = args[0]
    numSentences = len(os.listdir(minrule_grammars_loc))
    for iterNum in range(0, optsDict["iterations"]):
        #initialize counts to 0
        for lineNum in range(0, numSentences):
            minrule_fh = gzip.open(minrule_grammars_loc + "grammar.%d.gz"%(lineNum))
            sync_tree = tree(0, None, None, minrule_fh) #read in minimal grammar to tree struture
            #write function called updateCountsFromExample which takes in synctree
        print "hello"

if __name__ == "__main__":
    main()
