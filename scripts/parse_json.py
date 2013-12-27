#!/usr/bin/python -tt

'''
File: parse_json.py
Date: December 25, 2013
Description: This script looks at a JSON file that is
the output of the cdec decoder.  The JSON file describes
the hypergraph structure for each input sentence.  In 
particular, the hypergraph is a representation that describes
the intersection of the input sentence FSA and the SCFG,
to yield a sentence-specific grammar.  
Author: Avneesh Saluja (avneesh@cs.cmu.edu)
'''

import sys, commands, string, json, gzip
from pprint import pprint

def main():
    numSent = int(sys.argv[1])
    files_loc = sys.argv[2]
    for line_num in range(0, numSent):
        json_fh = gzip.open(files_loc + "%d.json.gz"%(line_num))
        data = json.load(json_fh)        
        print data.keys()
        node = data['node'] #node is a dict, where the keys are 'in_edges' and 'cat'
        in_edges = node['in_edges'] #in_edges is a list of IDs (rule or edge?)
        #print in_edges
        edges = data['edges']
        #for edge in edges:
        #    print edge
        pprint(data)
        json_fh.close()
        

if __name__ == "__main__":
    main()

