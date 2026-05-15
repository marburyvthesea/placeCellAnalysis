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
INPUT_numBins=${2:-32}

echo $INPUT_pathToAlignedData
echo "numBins=${INPUT_numBins}"

#add project directory to PATH
export PATH=$PATH/projects/p32501/


#load modules to use
module load matlab/r2023b

#cd to script directory
cd /home/ccc1839/placeCellAnalysis
#run analysis 

matlab -nosplash -nodesktop -r "addpath(genpath('/home/ccc1839/placeCellAnalysis'));maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));alignedFile='$INPUT_pathToAlignedData';numBins=str2double('$INPUT_numBins');run('placeCellAnalysisJJMquest.m');exit;"

echo 'finished analysis'
