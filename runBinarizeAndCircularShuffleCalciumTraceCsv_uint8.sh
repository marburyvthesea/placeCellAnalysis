#!/bin/bash
#SBATCH -A p30771
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -o ./logfiles/placeCellAnalysis.%x-%j.out # STDOUT
#SBATCH --job-name="binarizeShuffleCalciumUint8"
#SBATCH --mem=25G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16

module purge all

cd ~

# path to calcium trace CSV
INPUT_traceCsvPath=$1
INPUT_numShuffles=${2:-1000}
INPUT_binaryMarkMode=${3:-peak}
INPUT_outputDir=$4

echo $INPUT_traceCsvPath
echo "numShuffles=${INPUT_numShuffles}"
echo "binaryMarkMode=${INPUT_binaryMarkMode}"
echo "outputDir=${INPUT_outputDir}"

# add project directory to PATH
export PATH=$PATH/projects/p30771/

# load modules to use
module load matlab/r2023b

# cd to script directory
cd /home/jma819/placeCellAnalysis

# run analysis
matlab -nosplash -nodesktop -r "addpath(genpath('/home/jma819/placeCellAnalysis'));nCPUs=str2double(getenv('SLURM_CPUS_PER_TASK'));maxNumCompThreads(nCPUs);traceCsvPath='$INPUT_traceCsvPath';numShuffles=str2double('$INPUT_numShuffles');binaryMarkMode='$INPUT_binaryMarkMode';outputDir='$INPUT_outputDir';run('binarizeAndCircularShuffleCalciumTraceCsv_uint8.m');exit;"

echo 'finished analysis'
