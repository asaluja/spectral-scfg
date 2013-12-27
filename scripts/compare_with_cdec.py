#!/usr/bin/python -tt

'''
File: compare_with_cdec.py
Date: December 26, 2013
Description: this script takes the outputs of cdec's "show_target_graph" 
(edited by Avneesh to print out the relevant pieces of information)
and the output of intersect_scfg.py with the -d (debug) flag on, both of 
which contain information on the hypergraphs obtained after intersecting
the input sentence FSA with the SCFG, and looks at how they differ.  
Author: Avneesh Saluja (avneesh@cs.cmu.edu)
'''

import sys, commands, string

def readInRules(filename):
    file_list = open(filename, 'rb').readlines()
    line_num = 0
    nodeDict = {}
    while line_num < len(file_list):
        if file_list[line_num].split(':')[0] == "Node ID":            
            line_num += 1
            rules_used = []
            while (line_num < len(file_list)) and (file_list[line_num].split(':')[0] != "Node ID"): #until the next onee
                elements = file_list[line_num].strip().split(' ||| ')
                if len(elements) > 0:
                    rules_used.append(file_list[line_num].strip())
                line_num += 1
            rules_used = list(set(rules_used))
            #figure out span
            if len(rules_used) > 0:
                states = rules_used[0].strip().split(' ||| ')[0].split(',')
                span = states[1].split('-')
                cat = states[0]
                if cat == "[Goal":
                    continue
                key = (int(span[0]), int(span[1][:-1]))
                nodeDict[key] = rules_used
        else:
            line_num += 1
    return nodeDict        

def main():
    cdecDict = readInRules(sys.argv[1])
    ourDict = readInRules(sys.argv[2])
    for span in ourDict:
        if span not in cdecDict:
            print "Span (%d,%d) not in cdec rules!"%(span[0], span[1])
        else:
            setDiff = set(ourDict[span]) - set(cdecDict[span])
            if len(setDiff) > 0:
                print "Span (%d,%d) contain different rules for cdec output and our output"%(span[0], span[1])
                print "Differing rule(s) are: "
                for rule in setDiff:
                    print rule            
    for span in cdecDict:
        if span not in ourDict:
            print "Span (%d,%d) not in our rules!"%(span[0], span[1])
        else:
            setDiff = set(ourDict[span]) - set(cdecDict[span])
            if len(setDiff) > 0:
                print "Span (%d,%d) contain different rules for cdec output and our output"%(span[0], span[1])
                print "Differing rule(s) are: "
                for rule in setDiff:
                    print rule            
    print "Completed looking at differences"

if __name__ == "__main__":
    main()
