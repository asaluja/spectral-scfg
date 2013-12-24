#!/usr/bin/python -tt

import sys, commands, string

def main():
    feature_types = {}
    for line in sys.stdin:
        featureName = line.strip().split(':')[0]
        if len(featureName.split('_')) > 1:
            featureName = featureName.split('_')[0]
        if featureName not in feature_types:
            feature_types[featureName] = 1
        else:
            feature_types[featureName] += 1
    for feature in feature_types:
        print "%s:%d"%(feature, feature_types[feature])

if __name__ == "__main__":
    main()
