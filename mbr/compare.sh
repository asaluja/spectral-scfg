#!/bin/bash

x=$1

filt='cat'
if [ "$#" -eq 2 ]; then
  filt="grep $2"
fi

echo MLE-1
zcat /usr0/home/avneesh/spectral-scfg/data/evaluation/mle/mle-sentNorm-dev/grammar.$x.gz | ./convert-to-target-grammar.pl | ~/tools/cdec/training/utils/grammar_convert | ~/tools/cdec/decoder/cdec -f rescore -w weights.txt 2> /dev/null

echo
echo MLE-16
zcat /usr0/home/avneesh/spectral-scfg/data/evaluation/sentenceNorm/dec-rank16-dev/grammar.$x.gz | ./convert-to-target-grammar.pl  | ~/tools/cdec/training/utils/grammar_convert | ~/tools/cdec/decoder/cdec -f rescore -w weights.txt 2> /dev/null

echo
echo MLE-1
zcat /usr0/home/avneesh/spectral-scfg/data/evaluation/mle/mle-sentNorm-dev/grammar.$x.gz | ./convert-to-target-grammar.pl  | sort -k2 -t= -g | $filt | tail -15

echo
echo MLE-16
zcat /usr0/home/avneesh/spectral-scfg/data/evaluation/sentenceNorm/dec-rank16-dev/grammar.$x.gz | ./convert-to-target-grammar.pl  | sort -k2 -t= -g | $filt | tail -15
