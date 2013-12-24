#!/usr/bin/python -tt

import sys, commands, string

def main():
    for line in sys.stdin:
        line_format = line.strip().replace('(', '<')
        line_format = line_format.replace(')', '>')
        line_format = line_format.replace('@', '$')
        line_format = line_format.replace('[', '<')
        line_format = line_format.replace(']', '>')
        print line_format

if __name__ == "__main__":
    main()
