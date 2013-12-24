#!/usr/bin/python -tt
#########################
#Name: compare_sentence_marginals.py
#Date: September 27, 2013
#Description: given the output log files (stderr) of two systems' 
#parses on the same test data, we can compute how many times setup A 
#was better than setup B, and vice versa
#########################

import sys, commands, string

def populate_marginal_list(lines):
    marginals = []
    for i, line in enumerate(lines):
        if line == "Starting clock [parsing sentence]":
            prob_string_marg = lines[i+1].split("Parse weight: ")[1]
            if not prob_string_marg == "-Infinity":
                marginals.append(float(prob_string_marg))
    return marginals

def main():
    expA_lines = open(sys.argv[1], 'r').read().splitlines()
    expB_lines = open(sys.argv[2], 'r').read().splitlines()
    A_marginals = populate_marginal_list(expA_lines)
    B_marginals = populate_marginal_list(expB_lines)
    if len(A_marginals) != len(B_marginals):
        sys.stderr.write("Number of examples in file A (%d) and file B (%d) do not match!\n"%(len(A_marginals), len(B_marginals)))
        exit
    Abetter = 0
    Bbetter = 0
    same = 0
    for i, marginal in enumerate(A_marginals):
        if marginal > B_marginals[i]:
            Abetter += 1
        elif marginal < B_marginals[i]:
            Bbetter += 1
        else:
            same += 1
    print "Number of sentences where A > B: %d; B > A: %d; same: %d"%(Abetter, Bbetter, same)

if __name__ == "__main__":
    main()
