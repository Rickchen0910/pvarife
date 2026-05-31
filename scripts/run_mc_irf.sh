#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -l mem_free=16G
#$ -m be
#$ -M bc25911@essex.ac.uk
#$ -pe smp 40

BASE=/home/bc25911/MC_PVARIFE

mkdir -p $BASE/logs/success $BASE/logs/failed $BASE/csv $BASE/figures

/usr/bin/Rscript $BASE/mc_irf.R $1 $2 $3 $4
status=$?

echo "JOB_ID=$JOB_ID"
echo "Exit status=$status"

if [ $status -ne 0 ]; then
    echo "Job failed."
    mv $BASE/logs/out/${JOB_NAME}.o${JOB_ID} $BASE/logs/failed/ 2>/dev/null
    mv $BASE/logs/err/${JOB_NAME}.e${JOB_ID} $BASE/logs/failed/ 2>/dev/null
fi

exit $status
