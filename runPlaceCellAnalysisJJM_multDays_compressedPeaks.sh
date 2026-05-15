#!/bin/bash
#SBATCH -A p30771
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -o ./logfiles/placeCellAnalysis.%x-%j.out # STDOUT
#SBATCH --job-name="placeCellAnalysisMultiDaysCompressedPeaks"
#SBATCH --mem=25G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16

module purge all

cd ~

# path to combined GCAMP_with_velocity-style input file
INPUT_pathToAlignedData=$1
INPUT_numBins=${2:-32}

echo $INPUT_pathToAlignedData
echo "numBins=${INPUT_numBins}"

# add project directory to PATH
export PATH=$PATH/projects/p30771/

# load modules to use
module load matlab/r2023b

# cd to script directory
cd /home/jma819/placeCellAnalysis

# run analysis
matlab -nosplash -nodesktop -r "addpath(genpath('/home/jma819/placeCellAnalysis'));nCPUs=str2double(getenv('SLURM_CPUS_PER_TASK'));maxNumCompThreads(nCPUs);alignedFile='$INPUT_pathToAlignedData';numBins=str2double('$INPUT_numBins');run('placeCellAnalysisJJMquest_multDays_compressedPeaks.m');exit;"

echo 'finished analysis'
