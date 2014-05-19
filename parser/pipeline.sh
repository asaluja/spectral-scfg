#!/usr/bin/bash

scripts=$1
working=$2
params=$3
rank=$4
devSrc=$5
devSrcTgt=$6
devGrammar=$7
testSrc=$8
testSrcTgt=$9
testGrammar=${10}
config=${11}
numProc=${12}

#marginal computation, decoration of rules, for dev and devtest
#python ${scripts}/scripts/escape_special_characters.py < $devSrc > ${working}/dev.src.filt
cp $devSrc ${working}/dev.src
mkdir ${working}/rank${rank}-dev/
python ${scripts}/inside-outside/intersect_scfg.py -f $params $rank ${working}/dev.src $numProc ${working}/rank${rank}-dev/
mkdir ${working}/dec-rank${rank}-dev/
python ${scripts}/inside-outside/featurize_rules.py ${working}/rank${rank}-dev/ $devGrammar ${working}/dec-rank${rank}-dev/ $numProc
rm -rf ${working}/rank${rank}-dev/
python ${scripts}/scripts/corpus2sgm.py ${working}/dec-rank${rank}-dev/ < $devSrcTgt > ${working}/dev.sgm
rm ${working}/dev.src

#devtest
#python ${scripts}/scripts/escape_special_characters.py < $testSrc > ${working}/devtest.src.filt
cp $testSrc ${working}/devtest.src
mkdir ${working}/rank${rank}-devtest/
python ${scripts}/inside-outside/intersect_scfg.py -f $params $rank ${working}/devtest.src $numProc ${working}/rank${rank}-devtest/
mkdir ${working}/dec-rank${rank}-devtest/
python ${scripts}/inside-outside/featurize_rules.py ${working}/rank${rank}-devtest/ $testGrammar ${working}/dec-rank${rank}-devtest/ $numProc
rm -rf ${working}/rank${rank}-devtest/
python ${scripts}/scripts/corpus2sgm.py ${working}/dec-rank${rank}-devtest/ < $testSrcTgt > ${working}/devtest.sgm
rm ${working}/devtest.src

#MERT tuning
~/tools/cdec/training/dpmert/dpmert.pl --config $config --devset ${working}/dev.sgm --output-dir ${working}/rank${rank}.mert --weights ${working}/mert.weights.init --jobs $numProc

#evaluation
~/tools/cdec/training/utils/decode-and-evaluate.pl --jobs $numProc --input ${working}/dev.sgm --config $config --weights ${working}/rank${rank}.mert/weights.final &> rank${rank}.dev.bleu
~/tools/cdec/training/utils/decode-and-evaluate.pl --jobs $numProc --input ${working}/devtest.sgm --config $config --weights ${working}/rank${rank}.mert/weights.final &> rank${rank}.devtest.bleu
