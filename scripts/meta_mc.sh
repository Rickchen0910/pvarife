mkdir -p logs/out logs/err logs/success logs/failed csv

# Parameters to sweep
n_units=(25 50 100)
n_times=(25 50 100)
n_sds=500          # replications per (I,T) cell
n_cores=40         # must match -pe smp in run_mc_pvarife.sh

for nobs in ${n_units[@]}
do
  for tper in ${n_times[@]}
  do
    jobname=pvarife_mc_I${nobs}_T${tper}
    qsub -N $jobname \
         -o logs/out/${jobname}.o\$JOB_ID \
         -e logs/err/${jobname}.e\$JOB_ID \
         run_mc_pvarife.sh \
         $nobs $tper $n_sds
  done
done
