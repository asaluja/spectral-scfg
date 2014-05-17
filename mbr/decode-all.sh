#!/bin/bash

if [ "$#" -ne 1]; then
  echo Please give a path to grammar directory
  exit 1
fi

grammar_dir=$1
outf=output.$$
echo Decoding to $outf... 1>&2

#for x in `seq 0 2487`;
for x in `seq 0 60`;
do
  file=$grammar_dir/grammar.$x.gz
  echo !!!!!!!!!!!!! Processing $file 1>&2
  zcat $file | ./convert-to-target-grammar.pl | ~/tools/cdec/training/utils/grammar_convert | ~/tools/cdec/decoder/cdec -f rescore -w weights.txt >> $outf
done

