#!/usr/bin/bash -e

#######################
#File: decode_test.sh
#Date: December 27, 2013
#Author: Avneesh Saluja (avneesh@cs.cmu.edu)
#Description: for parallelized decoding of a test corpus, where we partition
#the test set into chunks and decode each chunk in a separate process.  The
#output of "intersect_scfg.py" is a set of all the rules in all of the parse
#trees/derivations of the source forest, where the NTs are decorated with the
#spans.  Each set of rules is written out to a sentence-specific file. 
######################

inputFile=$1
params=$2
rank=$3
numPartitions=$4
outDir=$5
script=/usr0/home/avneesh/spectral-scfg/code/inside-outside/intersect_scfg.py

for (( i = 0; i < numPartitions; i++ )); do
    python $script $params $inputFile $numPartitions $i $rank $outDir &
done
wait
echo "completed marginal computations for $inputFile"
    