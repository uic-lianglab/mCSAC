
#$ -N J.q.463113
#$ -S /bin/bash
#$ -o Job.log
#$ -j y
#$ -m n
#$ -l h_rt=00:00:50
#$ -cwd
../build/bin/ensemble -conf ./ensemble_try.ini -prefix _H -samplesize 5

