config = struct();
config.storageMode = 'compressed_uint8';
config.outputSuffix = '_placeCellAnalysis_multDays_compressedPeaks';
config.peakDeflateLevel = 9;
config.peakChunkShuffles = 1;
if exist('numBins', 'var') && ~isempty(numBins)
    config.numBins = numBins;
end
if exist('runTimestamp', 'var') && ~isempty(runTimestamp)
    config.runTimestamp = runTimestamp;
end

placeCellAnalysisJJMquest_multDays_savedShuffles(alignedFile, config);
