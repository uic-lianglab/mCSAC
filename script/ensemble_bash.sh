#!/bin/bash - 

Prefix=($(echo {A..Z}{A..Z}{A..Z}))
#Prefix=($(echo {K..Z}{A..Z}))

# Ensemble size = 100000
SampleSize=1
idx=0
while [ $idx -lt 100 ]
do
	../build/bin/ensemble -conf ensemble_try.ini -prefix ${Prefix[$idx]} -samplesize $SampleSize

idx=$(( $idx + 1 ))
done

