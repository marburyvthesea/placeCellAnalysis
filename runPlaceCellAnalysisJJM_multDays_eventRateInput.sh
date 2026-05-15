#!/bin/bash
#SBATCH -A p30771
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -o ./logfiles/placeCellAnalysis.%x-%j.out # STDOUT
#SBATCH --job-name="placeCellAnalysisMultiDaysEventRateInput"
#SBATCH --mem=25G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16

module purge all

cd ~

# path to combined instantaneousEventRate-style input file
INPUT_pathToEventRateData=$1
INPUT_numBins=${2:-32}
INPUT_pathToMetadataData=$3

echo $INPUT_pathToEventRateData
echo "numBins=${INPUT_numBins}"
echo "metadataFile=${INPUT_pathToMetadataData}"

# add project directory to PATH
export PATH=$PATH/projects/p30771/

# load modules to use
module load matlab/r2023b

# cd to script directory
cd /home/jma819/placeCellAnalysis

# run analysis
matlab -nosplash -nodesktop -r "addpath(genpath('/home/jma819/placeCellAnalysis'));nCPUs=str2double(getenv('SLURM_CPUS_PER_TASK'));maxNumCompThreads(nCPUs);alignedFile='$INPUT_pathToEventRateData';numBins=str2double('$INPUT_numBins');metadataFile='$INPUT_pathToMetadataData';run('placeCellAnalysisJJMquest_multDays_eventRateInput.m');exit;"

echo 'finished analysis'
