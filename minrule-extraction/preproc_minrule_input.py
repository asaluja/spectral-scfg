#!/usr/bin/python -tt

'''
File: preproc_minrule_input.py
Date: September 30, 2013
Description: this script takes an alignment file, in the standard format, e.g., 0-0 1-0 1-1 2-2 etc.
and changes the format to what is required for the minimal rule extraction code.
Author: Avneesh Saluja (avneesh@cs.cmu.edu)
'''

import sys, commands, string

def convertTupleListToDict(list_of_tuples):
    dictToReturn = {}
    for item in list_of_tuples:
        if item[0] not in dictToReturn: 
            index_list = []
            index_list.append(item[1])
            dictToReturn[item[0]] = index_list
        else:
            dictToReturn[item[0]].append(item[1])
    return dictToReturn

def main():
    for line in sys.stdin:
        alignments = [(int(alignment.split('-')[0]), int(alignment.split('-')[1])) for alignment in line.strip().split()]
        numSource = max([alignment[0] for alignment in alignments]) + 1
        numTarget = max([alignment[1] for alignment in alignments]) + 1
        print "%d %d"%(numSource, numTarget) #print dimensions of alignment matrix
        alignDict = convertTupleListToDict(alignments)
        for srcIdx in range(0, numSource):
            if srcIdx not in alignDict:
                print
            else:
                print ' '.join([str(x+1) for x in alignDict[srcIdx]])

if __name__ == "__main__":
    main()
