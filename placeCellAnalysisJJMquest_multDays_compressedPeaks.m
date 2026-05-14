config = struct();
config.storageMode = 'compressed_uint8';
config.outputSuffix = '_placeCellAnalysis_multDays_compressedPeaks';
config.peakDeflateLevel = 9;
config.peakChunkShuffles = 1;

placeCellAnalysisJJMquest_multDays_savedShuffles(alignedFile, config);
