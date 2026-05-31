#!/bin/bash
#SBATCH --job-name=pvarife_mc
#SBATCH --array=1-9                 # 9 tasks: I in {25,50,100} x T in {25,50,100}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16          # cores per task (MonteCarlo ncpus)
#SBATCH --mem=8G
#SBATCH --time=04:00:00             # 4h per task (generous for n_sds=500, n_out=50)
#SBATCH --output=logs/mc_%A_%a.out  # stdout: logs/mc_<jobid>_<taskid>.out
#SBATCH --error=logs/mc_%A_%a.err   # stderr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bc25911@essex.ac.uk

# ---------------------------------------------------------------------------
# Essex VERNE / your HPC: uncomment and set the correct R module
# module purge
# module load R/4.3.1
# ---------------------------------------------------------------------------

mkdir -p logs mc_hpc_results

echo "=== Task $SLURM_ARRAY_TASK_ID started at $(date) ==="
echo "Node: $SLURMD_NODENAME  |  Cores: $SLURM_CPUS_PER_TASK"

Rscript mc_convergence_hpc.R \
  --task   "$SLURM_ARRAY_TASK_ID" \
  --ncores "$SLURM_CPUS_PER_TASK" \
  --nsds   500 \
  --outdir mc_hpc_results

echo "=== Task $SLURM_ARRAY_TASK_ID finished at $(date) ==="
