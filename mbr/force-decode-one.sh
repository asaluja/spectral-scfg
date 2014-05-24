#!/bin/bash

scriptLoc=/usr0/home/avneesh/spectral-scfg/code/mbr
working=/usr0/home/avneesh/spectral-scfg/data/evaluation/sentenceNorm/force-decode
grammarFile=$1
echo $grammarFile
parallel=/usr0/home/avneesh/spectral-scfg/data/corpus/training/training.filter.no-invalid.de-en
x=`echo $grammarFile | rev | cut -d'.' -f2`
y=$(($x + 1))

zcat $grammarFile | ${scriptLoc}/convert-to-target-grammar.pl | ~/tools/cdec/training/utils/grammar_convert > ${working}/${x}.hg
head -$y  $parallel | tail -1 | ~/tools/cdec/corpus/cut-corpus.pl 2 > ${working}/ref.${x}

paste ${working}/${x}.hg ${working}/ref.${x} | perl -e 'while(<>){s/\t/ ||| /;print;}' | ~/tools/cdec/decoder/cdec -z -f rescore -w ${working}/weights.txt
rm ${working}/${x}.hg
rm ${working}/ref.${x}

