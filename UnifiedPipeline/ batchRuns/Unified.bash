#!/bin/bash
##Defines Job Array Properties
#SBATCH --job-name=myshortrun
#SBATCH --time=08:00:00
#SBATCH --partition=short
#SBATCH --ntasks=10
#SBATCH --nodes=1
###SBATCH --gres=gpu:1
#SBATCH --mem=32Gb
#SBATCH --output=GN%a.out
##SBATCH --error=GN%a.err
#SBATCH --array=1-1000%1000

##Loads dependencies
module load cuda/11.3
module load anaconda3
source activate blind

##Reads User Input
echo "TotJobs: $1"
echo "jobStart: $2"
echo "inputFile: $3"

##Defines environmental variables
stepNum=$(($SLURM_ARRAY_TASK_ID + $2))
parentdir="$(dirname "$PWD")"

# Runs simulation
$HOME/.conda/envs/AIBind/bin/python3.9 $parentdir/blind.py --instance $stepNum --total_jobs $1 --file $3
