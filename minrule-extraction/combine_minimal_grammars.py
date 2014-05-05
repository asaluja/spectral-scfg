#!/usr/bin/python -tt

'''
File: combine_minimal_grammars.py
Date: April 29, 2014
Description: this script takes as input the list of invalid
rules (from the STDERR output of tree_to_rule.py), as well 
as the input directory (list of .gz rules produced by the
tree_to_rule.py script) and either copies over the original
.gz file produced by the tree_to_rule script, or if the sentence
is invalid, uses the other minimal grammar sentence. 
Arg1: list of invalid rules
Arg2: directory where guaranteed binary minimal rules are
Arg3: directory where ZGC rules are
Arg4: output directory
Usage: 
'''

import sys, commands, string, os, gzip

def ProcessInvalidRules(filename):
    invalid_rules = {}
    for line in open(filename, 'rb'):
        sentNum_Str = line.strip().split(':')[0]
        sentNum = int(sentNum_Str.split()[1])
        invalid_rules[sentNum] = 1
    return invalid_rules
    

def main():
    invalid_rules = ProcessInvalidRules(sys.argv[1])
    binary_grammar_loc = sys.argv[2]
    minrule_grammar_loc = sys.argv[3]
    out_grammar_loc = sys.argv[4]
    numFiles = len(os.listdir(minrule_grammar_loc))
    for x in xrange(0, numFiles):
        file_to_copy = binary_grammar_loc + "/grammar.%d.gz"%x if x in invalid_rules else minrule_grammar_loc + "/grammar.%d.gz"%x
        out_filename = out_grammar_loc + "/grammar.%d.gz"%x
        fh_in = gzip.open(file_to_copy, 'rb')
        fh_out = gzip.open(out_filename, 'w')
        for rule in fh_in:
            fh_out.write("%s"%rule)
        fh_out.close()
        fh_in.close()

if __name__ == "__main__":
    main()
