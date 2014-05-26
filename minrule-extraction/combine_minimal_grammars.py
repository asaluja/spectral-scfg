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

import sys, commands, string, os, gzip, re

def ProcessInvalidRules(filename):
    invalid_rules = {}
    for line in open(filename, 'rb'):
        sentNum_Str = line.strip().split(':')[0]
        sentNum = int(sentNum_Str.split()[1])
        invalid_rules[sentNum] = 1
    return invalid_rules

def ReformatRule(rule):
    elements = rule.split(' ||| ')
    expr = re.compile(r'\[([^]]+)\]')
    new_src = []
    NT_counter = 1
    for item in elements[1].split():
        if expr.match(item):
            new_src.append("[X,%d]"%NT_counter)
            NT_counter += 1
        else:
            new_src.append(item)
    new_src_str = ' '.join(new_src)
    new_tgt = []
    for item in elements[2].split():
        if expr.match(item):
            NT_counter = int(re.findall(expr, item)[0])
            new_tgt.append("[X,%d]"%NT_counter)
        else:
            new_tgt.append(item)
    new_tgt_str = ' '.join(new_tgt)
    return ' ||| '.join([elements[0], new_src_str, new_tgt_str])

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
            rule_format = ReformatRule(rule.strip()) if x in invalid_rules else rule.strip()
            fh_out.write("%s\n"%rule_format)
        fh_out.close()
        fh_in.close()

if __name__ == "__main__":
    main()
