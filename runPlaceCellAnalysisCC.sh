#!/bin/bash
#SBATCH -A p32501
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH -o ./logfiles/placeCellAnalysis.%x-%j.out # STDOUT
#SBATCH --job-name="placeCellAnalysis"
#SBATCH --mem-per-cpu=5200M
#SBATCH -N 1
#SBATCH -n 16

module purge all

cd ~

#path to file 

INPUT_pathToAlignedData=$1

echo $INPUT_pathToAlignedData

#add project directory to PATH
export PATH=$PATH/projects/p32501/


#load modules to use
module load matlab/r2023b

#cd to script directory
cd /home/ccc1839/placeCellAnalysis
#run analysis 

matlab -nosplash -nodesktop -r "addpath(genpath('/home/ccc1839/placeCellAnalysis'));maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));alignedFile='$INPUT_pathToAlignedData';run('placeCellAnalysisJJMquest.m');exit;"

echo 'finished analysis'
