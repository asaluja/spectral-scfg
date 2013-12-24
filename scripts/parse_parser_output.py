#!/usr/bin/python -tt

#########################
#File: parse_parser_output.py
#Date: September 2013
#Description: parses the stderr/log of a parse
#on test sentences, and computes log likelihoods, 
#log marginals, and number of examples parsed
#########################

import sys, commands, string, math

def main():
    lines = sys.stdin.read().splitlines()
    logprobs = 0
    marginals = 0
    numEntries = 0
    numTotal = 1
    for i,line in enumerate(lines):
        if line == "Starting clock [parsing sentence]":
            prob_string_marg = lines[i+1].split("Parse weight: ")[1]
            if not prob_string_marg == "-Infinity":
                marginals += math.log(float(prob_string_marg))
            else:
                sys.stderr.write("%d\n"%(numTotal))
            prob_string_io = lines[i+2].split(": +[")[-1].rstrip(']')
            if not prob_string_io == "-Infinity":
                logprobs += math.log(float(prob_string_io))
                numEntries += 1
            numTotal += 1
    print "Log Likelihood: %.3f; Log Marginal: %.3f; Number of examples: %d"%(logprobs, marginals, numEntries)

if __name__ == "__main__":
    main()
