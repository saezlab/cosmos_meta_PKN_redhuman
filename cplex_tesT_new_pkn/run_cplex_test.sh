#!/usr/local_rwth/bin/zsh
#
#SBATCH --job-name=cplex_test
#SBATCH --output=cplex_test.out
#SBATCH --error=cplex_test.err
#SBATCH --export=ALL
#SBATCH --gpus=0
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --no-requeue
#SBATCH --mem=180G
#SBATCH -t 30:00:00
#SBATCH --account=jrc_combine

Rscript cplex_test.R