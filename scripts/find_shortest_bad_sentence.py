#!/usr/bin/python -tt

import sys, commands, string

bad_sentences = map(int, open(sys.argv[1], 'r').read().splitlines())
line_counter = 1
shortest_sentence_length = 1000
shortest_sentence_number = -1
for line in sys.stdin:
    if line_counter in bad_sentences:
        sent_length = len(line.strip().split(' | '))
        if sent_length < shortest_sentence_length and sent_length > 4:
            shortest_sentence_number = line_counter
            shortest_sentence_length = sent_length
    line_counter += 1
print "Shortest sentence is sentence number %d"%(shortest_sentence_number)
