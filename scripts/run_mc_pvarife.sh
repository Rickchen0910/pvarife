#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -l mem_free=8G
#$ -m be
#$ -M bc25911@essex.ac.uk
#$ -pe smp 40

BASE=/home/bc25911/MC_PVARIFE

mkdir -p $BASE/logs/success
mkdir -p $BASE/logs/failed
mkdir -p $BASE/csv

/usr/bin/Rscript $BASE/mc_pvarife.R $1 $2 $3
status=$?

echo "JOB_ID=$JOB_ID"
echo "Exit status=$status"

if [ $status -ne 0 ]; then
    echo "Job failed. Moving logs..."
    mv $BASE/logs/out/${JOB_NAME}.o${JOB_ID} \
       $BASE/logs/failed/ 2>/dev/null
    mv $BASE/logs/err/${JOB_NAME}.e${JOB_ID} \
       $BASE/logs/failed/ 2>/dev/null
fi

exit $status
