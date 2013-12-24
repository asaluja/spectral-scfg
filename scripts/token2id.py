#!/usr/bin/python -tt

import sys, commands, string

def main():
    token2ID = {}
    tokenID = 0
    delimiter = ""
    if len(sys.argv) > 2:
        delimiter = sys.argv[2] #e.g.: '|'; default is space
    else:
        delimiter = " "
    outF = open(sys.argv[1], 'w')
    for line in sys.stdin:
        line_format = []
        tokens = line.strip().split(delimiter)
        for token in tokens: #token may be a phrase, word, POS sequence, or POS
            token = token.rstrip().lstrip()
            if token not in token2ID: #phrase doesn't exist
                token2ID[token] = tokenID
                print >> outF, "%d:%s"%(tokenID, token)
                tokenID += 1
            line_format.append(token2ID[token])
        print '|'.join(map(str, line_format))
    outF.close()
    print >> sys.stderr, "Number of unique units: %d"%(len(token2ID))
    maxL = 0
    longestToken = ""
    for token in token2ID.keys():
        if len(token.split()) > maxL:
            maxL = len(token.split())
            longestToken = token
    print >> sys.stderr, "Length of longest unit is: %d"%(maxL)
    print >> sys.stderr, "Unit is: %s"%(longestToken)


if __name__ == "__main__":
    main()
