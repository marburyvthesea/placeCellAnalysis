config = struct();
config.storageMode = 'sparse_linear_idx';
config.outputSuffix = '_placeCellAnalysis_multDays_sparseIndices';
config.peakDeflateLevel = 9;
config.indexChunkEntries = 1000000;
if exist('numBins', 'var') && ~isempty(numBins)
    config.numBins = numBins;
end
if exist('runTimestamp', 'var') && ~isempty(runTimestamp)
    config.runTimestamp = runTimestamp;
end

placeCellAnalysisJJMquest_multDays_savedShuffles(alignedFile, config);
