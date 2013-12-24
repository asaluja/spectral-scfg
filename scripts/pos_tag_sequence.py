#!/usr/bin/python -tt

############################
#File: pos_tag_sequence.py
#Date: September 12, 2013
#Description: code using NLTK to POS-tag a sequence
#The code assumes a particular input format (phrase1 | phrase2 | phrase3 etc)
#It strips the delimiters, tags the sentence, and writes out the tags with
#the delimiters put back in.
#Author: Avneesh Saluja (avneesh@cs.cmu.edu)
###########################

import sys, commands, string, nltk

def main():
    for line in sys.stdin:
        line_no_delim_list = ' '.join(line.strip().split('|')).split() #for the specific format of our data, remove the | and replace with space
        tag_sequence = nltk.pos_tag(line_no_delim_list) #tag the words
        tag_sequence.append(('|', '|')) #add the delimiter back in
        tag_dict = {}
        #N.B: code below will fail if we have one type, two tokens in a sentence and the tokens have different POS
        for k, v in tag_sequence: #convert list to dict
            tag_dict[k] = v
        pos_tag_line = []
        for token in line.split(): #go through all tokens in sentence, including delimiter '|'
            pos_tag_line.append(tag_dict[token]) #convert to POS tag
        print ' '.join(pos_tag_line) #print POS sequence out as string

if __name__ == "__main__":
    main()
