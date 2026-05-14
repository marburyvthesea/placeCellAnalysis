config = struct();
config.storageMode = 'sparse_linear_idx';
config.outputSuffix = '_placeCellAnalysis_multDays_sparseIndices';
config.peakDeflateLevel = 9;
config.indexChunkEntries = 1000000;

placeCellAnalysisJJMquest_multDays_savedShuffles(alignedFile, config);
