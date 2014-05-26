#!/bin/bash

scriptLoc=/usr0/home/avneesh/spectral-scfg/code/mbr
grammarLoc=$1
parallel=$2
working=$3
outputFile=$4
rm $outputFile #if already exists, remove it
max_jobs=$5
echo "max jobs: $max_jobs"

numFiles=`ls -1 $grammarLoc | wc -l`
echo "Number of files: $numFiles"
seqLimit=$(($numFiles-1))

seq 0 $seqLimit | xargs  -i --max-procs=$max_jobs bash -c "${scriptLoc}/force-decode-slave.sh $working $parallel $grammarLoc {}; echo {} done"
wait
for (( i = 0; i < $numFiles; i++ )); do
    norm=`grep 'Partition         log(Z)' ${working}/${i}.ll | cut -d':' -f2 | sed 's/^[[:space:]]*//'`
    unreachable=`grep 'REFERENCE UNREACHABLE' ${working}/${i}.ll | wc -l`
    if [[ $unreachable -eq 1 ]]; then
	echo "LL: 0|$norm" >> $outputFile
    else
	#logp=`grep 'Constr. forest  Viterbi logp' ${working}/${i}.ll | cut -d':' -f2 | sed 's/^[[:space:]]*//'`
	logp=`grep 'Contst. partition  log(Z)'  ${working}/${i}.ll | cut -d':' -f2 | sed 's/^[[:space:]]*//'`	
	echo "LL: $logp|$norm" >> $outputFile
    fi
    rm ${working}/${i}.ll
done