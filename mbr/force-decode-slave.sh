#!/bin/bash

scriptLoc=/usr0/home/avneesh/spectral-scfg/code/mbr
working=$1
parallel=$2
grammarLoc=$3
#grammarFile=$3

#x=`echo $grammarFile | rev | cut -d'.' -f2`
x=$4
y=$(($x + 1))

zcat ${grammarLoc}/grammar.$x.gz | ${scriptLoc}/convert-to-target-grammar.pl | ~/tools/cdec/training/utils/grammar_convert > ${working}/${x}.hg
head -$y  $parallel | tail -1 | ~/tools/cdec/corpus/cut-corpus.pl 2 > ${working}/ref.${x}

paste ${working}/${x}.hg ${working}/ref.${x} | perl -e 'while(<>){s/\t/ ||| /;print;}' | ~/tools/cdec/decoder/cdec -z -f rescore -w ${working}/weights.txt &> ${working}/${x}.ll
rm ${working}/${x}.hg
rm ${working}/ref.${x}