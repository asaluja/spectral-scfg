import sys, os, getopt

def treefy_tags(line):
    tokens = line.split('|')
    head = '('
    tail = ')'

    for token in tokens:
        if tagdef:
            (tag, word) = token.rsplit('/', 1)
        else:
            (phrase, tag) = (token, "1")
        phrase = phrase.rstrip().lstrip()
        head += '(' + tag + ' '
        head += '(' + tag + 'LEFT ' #not sure if we need this - if we don't, check about closing brackets
        #head += '(' + raw + 'LEFT%' + position + ' '
        head += phrase + ')' + ' '
        tail += ')'
    #head += '(STOP%!!!E STOP)'
    head += '(STOP STOP)'
    return head + tail

def tagfy_tree(line):
    tokens = line.split()
    seq = ''
    for i, token in enumerate(tokens[:-1]):
        if ')' in token:
            (word, junk) = token.split(')', 1)
            tag = tokens[i-1][1:-4]
            if tagdef:
                seq += tag + '/' + word + ' '
            else:
                seq += word + '/' + tag + ' '
    return seq 

seq2tree = True
tagdef = False
(opts, args) = getopt.getopt(sys.argv[1:], 'ro')
for opt in opts:
    if opt[0] == '-r':
        seq2tree = False
    if opt[0] == '-o':
        tagdef = True
infile = open(args[0])

lines = infile.readlines()
for line in lines:
    if line != '\n':
        if seq2tree:
            print treefy_tags(line)
        else:
            print tagfy_tree(line)
