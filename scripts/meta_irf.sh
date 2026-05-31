mkdir -p logs/out logs/err logs/success logs/failed csv figures

# Parameters — match the paper: I in {25,50,100}, T in {25,50,100}
n_units=(25 50 100)
n_times=(25 50 100)
n_sds=500
# identification: "sr" = short-run (Figure S.1), "lr" = long-run (Figure S.2)
id="sr"

for nobs in ${n_units[@]}
do
  for tper in ${n_times[@]}
  do
    jobname=pvarife_irf_${id}_I${nobs}_T${tper}
    qsub -N $jobname \
         -o logs/out/${jobname}.o\$JOB_ID \
         -e logs/err/${jobname}.e\$JOB_ID \
         run_mc_irf.sh \
         $nobs $tper $n_sds $id
  done
done
