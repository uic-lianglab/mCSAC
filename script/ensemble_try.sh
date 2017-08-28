#!/bin/bash - 

Prefix=($(echo {R..Z}{A..Z}{A..Z}{A..Z}))
#Prefix2=($(echo {A..Z}{A..Z}))

# Ensemble size = 9000
SampleSize=5
idx=0
while [ $idx -lt 50000 ]
do
  ItemID=${Prefix[$idx]}
  FJOB=$(basename $ItemID).q
  JOBNAME=J$FJOB.$idx
  FileLog=Job.log

echo "
#$ -N $JOBNAME
#$ -S /bin/bash
#$ -o $FileLog
#$ -j y
#$ -m n
#$ -l h_rt=00:00:50
#$ -cwd
../build/bin/ensemble -conf ./ensemble_try.ini -prefix ${Prefix[$idx]}_A -samplesize $SampleSize
" >  $FJOB
qsub $FJOB
rm -f $FJOB
idx=$(( $idx + 1 ))
done

