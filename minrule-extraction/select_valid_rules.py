#!/usr/bin/python -tt

'''
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
    minrule_grammar_loc = sys.argv[2]
    out_grammar_loc = sys.argv[3]
    numFiles = len(os.listdir(minrule_grammar_loc))
    counter = 0
    for x in xrange(0, numFiles):
        if x not in invalid_rules:            
            file_to_copy = minrule_grammar_loc + "/grammar.%d.gz"%x
            out_filename = out_grammar_loc + "/grammar.%d.gz"%counter
            fh_in = gzip.open(file_to_copy, 'rb')
            fh_out = gzip.open(out_filename, 'w')
            counter += 1
            for rule in fh_in:
                fh_out.write("%s"%rule)
            fh_out.close()
            fh_in.close()
if __name__ == "__main__":
    main()
