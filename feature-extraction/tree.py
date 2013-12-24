#!/usr/bin/python -tt

import sys, commands, string, gzip, re

class tree:
    '''
    Tree class is a recursive data structure used for feature extraction.  
    Since the features we use are based on the tree structure, it's 
    important to convert the per sentence minimal grammar to the tree 
    so that we can extract these features easier in the main function 
    of the module. 
    '''
    def __init__(self, root_idx, root, rule_dict, filehandle=None):
        if filehandle:
            rule_dict = self.read_grammar_to_dict(filehandle)
        raw_rule = rule_dict[root_idx]
        pattern = re.compile('\[[0-9]+')        
        self.rule = pattern.sub('[X', raw_rule)
        if root == None: #i.e., we are at the root of the sentence
            self.rule = self.rule.replace('[X]', '[S]')
        self.src = self.rule.split(' ||| ')[1] #extract src and tgt phrases
        self.tgt = self.rule.split(' ||| ')[2]
        self.children = []
        self.parent = root
        children_idx = [int(nt.split(',')[0]) for nt in re.findall(r'\[([^]]+)\]', raw_rule.split(' ||| ')[1])] #extract NTs to find children
        for child in children_idx: #recurse on children
            child_tree = tree(child, self, rule_dict)
            self.children.append(child_tree) #once the child structure has been set up, add to current node's children

    '''
    The per-sentence grammar is stored in a dict, with each NT given a unique
    symbol, and the actual rule as the value. 
    '''
    def read_grammar_to_dict(self, filehandle):
        rule_dict = {}
        for rule in filehandle:
            k = rule.strip().split(' ||| ')[0]
            raw_symbol = int(re.findall(r'\[([^]]+)\]', k)[0]) #raw symbol is a number
            rule_dict[raw_symbol] = rule.strip()
        return rule_dict


            


            
