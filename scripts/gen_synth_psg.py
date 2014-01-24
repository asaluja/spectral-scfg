#!/usr/bin/python -tt

'''
File: gen_synth_psg.poy
Date:
Description: 
'''

import sys, commands, string, os, gzip

def main():
    dirOut = sys.argv[1]
    if not os.path.exists(dirOut):
        os.makedirs(dirOut)        
    for count,line in enumerate(sys.stdin):
        sentFile = gzip.open(dirOut + "/grammar.%d.gz"%(count), 'wb')
        elements = line.strip().split(' ||| ')
        if elements[0] == 'b':
            sentFile.write("[0] ||| [1,1] ||| [1,1]\n")
            sentFile.write("[1] ||| b ||| BB\n")
        else:
            sentFile.write("[0] ||| [1,1] [2,2] ||| [1,1] [2,2]\n")
            sentFile.write("[1] ||| a ||| A\n")
            sentFile.write("[2] ||| b ||| B\n")
        sentFile.close()

if __name__ == "__main__":
    main()
