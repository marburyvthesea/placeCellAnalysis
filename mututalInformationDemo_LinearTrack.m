% spatial_mutual_information_linear_track_demo.m
%
% Synthetic demonstration of why occupancy-weighted spatial mutual
% information is safer than using only a high conditional firing
% probability P(event | bin) to classify place cells on a nearly
% 2-dimensional linear track.
%
% The "track" here is an elongated corridor with a small but nonzero
% vertical width, so the geometry is still 2D while behaving mostly like a
% linear track.
%
% Metric implemented here:
%   I_pos(x_i) = sum_{k in {0,1}} P(k | x_i) * log2( P(k | x_i) / P(k) )
%   SI         = sum_i P(x_i) * I_pos(x_i)
%
% This is the binary-event spatial mutual information used in:
%   Kinsky et al., Current Biology (2018)
%   "Hippocampal Place Fields Maintain a Coherent and Flexible Map across
%    Long Timescales", doi: 10.1016/j.cub.2018.09.037
%
% What this script shows:
%   1) A true place field on the track.
%   2) A fast pass-through hotspot: the animal moves quickly through one
%      track segment, so occupancy there is tiny. Peak P(event | bin) looks
%      strong, but SI stays small because P(x) is tiny.
%   3) A slow linger zone: the animal spends a long time in one zone, but a
%      mostly uniform high-baseline cell only has a tiny local bump up to
%      the same peak P(event | bin). Peak P(event | bin) is high, yet SI
%      stays small because almost all bins look similar to the global rate.
%
% Output:
%   - spatial_mi_linear_track_maps.png
%   - spatial_mi_linear_track_summary.png
%
% Notes:
%   - This is a pedagogical simulation, not a re-analysis of the paper.
%   - The place-cell significance test uses circularly shifted events.

clear;
clc;
close all;

rng(11);

%% User-tunable parameters
nx = 72;                    % bins along the long axis of the track
ny = 16;                    % bins across the width of the track
nFrames = 4000;             % synthetic imaging frames
nShuffles = 1000;           % paper-style shuffle count
naivePeakThreshold = 0.50;  % intentionally naive "high P(event|bin)" rule
peakTarget = 0.60;          % all synthetic cells share this peak probability
outDir = pwd;

%% Build a narrow 2D linear track
[xGrid, yGrid] = meshgrid(1:nx, 1:ny);
trackMidY = (ny + 1) / 2;

% A slight waviness keeps the geometry "nearly 2D" rather than perfectly flat.
trackCenterline = trackMidY + 0.75 * sin(2 * pi * (xGrid - 1) / (nx - 1));
trackHalfWidth = 1.8 + 0.25 * cos(2 * pi * (xGrid - 1) / (nx - 1) + 0.35);
trackMask = abs(yGrid - trackCenterline) <= trackHalfWidth;

trackProfile = exp(-((yGrid - trackCenterline) .^ 2) ./ (2 * 0.95 ^ 2)) .* trackMask;
baseLongitudinal = 0.88 + 0.12 * exp(-((xGrid - 0.52 * nx) .^ 2) ./ (2 * 16 ^ 2));
baseOccMap = trackProfile .* baseLongitudinal;
baseOccMap = normalizeProbMap(baseOccMap);

%% Velocity/occupancy edge cases
fastTransitX = 15;
slowLingerX = 58;
trueFieldX = 30;

% Increase 0.97 to make the fast-transit segment even less occupied.
fastTransitOccMap = baseOccMap .* ...
    (1 - 0.97 * maxNormalize(trackAlignedField(xGrid, yGrid, trackCenterline, fastTransitX, 2.8, 0.9))) ...
    + 1e-7 * trackMask;
fastTransitOccMap = normalizeProbMap(fastTransitOccMap);

% Increase 1.9 to make the slow-linger zone even more occupied.
% Keep the added occupancy on the track itself; otherwise off-track bins
% acquire occupancy but still have P(event|x)=0, which artificially
% inflates spatial information.
slowLingerOccMap = baseOccMap + ...
    1.9 * trackAlignedField(xGrid, yGrid, trackCenterline, slowLingerX, 6.0, 1.2) .* trackMask;
slowLingerOccMap = normalizeProbMap(slowLingerOccMap);

%% Synthetic conditional event-probability maps
placeMap = 0.02 + (peakTarget - 0.02) * ...
    maxNormalize(trackAlignedField(xGrid, yGrid, trackCenterline, trueFieldX, 4.2, 0.9));
placeMap = placeMap .* trackMask;

fastTransitHotspotMap = 0.02 + (peakTarget - 0.02) * ...
    maxNormalize(trackAlignedField(xGrid, yGrid, trackCenterline, fastTransitX, 1.4, 0.7));
fastTransitHotspotMap = fastTransitHotspotMap .* trackMask;

% Keep the peak at 0.60, but make the modulation tiny and very local on
% top of a nearly uniform high baseline. This creates a cleaner "naive peak
% says yes, SI says no" edge case on the linear track.
slowLingerBaseline = 0.52;
slowLingerStateMap = slowLingerBaseline + (peakTarget - slowLingerBaseline) * ...
    maxNormalize(trackAlignedField(xGrid, yGrid, trackCenterline, slowLingerX, 1.3, 0.55));
slowLingerStateMap = slowLingerStateMap .* trackMask;

scenarios = struct([]);
scenarios(1).name = 'True track field';
scenarios(1).shortName = 'True field';
scenarios(1).occMap = baseOccMap;
scenarios(1).pEventGivenX = placeMap;
scenarios(1).description = 'Compact field on a narrow track';

scenarios(2).name = 'Fast pass-through segment';
scenarios(2).shortName = 'Fast transit';
scenarios(2).occMap = fastTransitOccMap;
scenarios(2).pEventGivenX = fastTransitHotspotMap;
scenarios(2).description = 'High speed produces a low-occupancy hotspot';

scenarios(3).name = 'Slow linger zone';
scenarios(3).shortName = 'Slow linger';
scenarios(3).occMap = slowLingerOccMap;
scenarios(3).pEventGivenX = slowLingerStateMap;
scenarios(3).description = 'Long dwell time plus only a tiny bump on a high baseline';

nScenarios = numel(scenarios);

%% Expected SI from the underlying maps
for i = 1:nScenarios
    [expectedSI, localInfoMap, pEventGlobal] = spatialMIFromMaps( ...
        scenarios(i).occMap, scenarios(i).pEventGivenX);

    [~, peakIdx] = max(scenarios(i).pEventGivenX(:));

    scenarios(i).expectedSI = expectedSI;
    scenarios(i).localInfoMap = localInfoMap;
    scenarios(i).globalEventProb = pEventGlobal;
    scenarios(i).peakConditional = max(scenarios(i).pEventGivenX(:));
    scenarios(i).peakBinOccupancy = scenarios(i).occMap(peakIdx);
    scenarios(i).trackOccupancyMarginal = sum(scenarios(i).occMap, 1);
    scenarios(i).trackEventMarginal = weightedTrackAverage( ...
        scenarios(i).pEventGivenX, scenarios(i).occMap);
    scenarios(i).trackLocalSIMarginal = sum(scenarios(i).localInfoMap, 1);
end

%% Simulate one session per cell and evaluate the paper-style shuffle rule
for i = 1:nScenarios
    [posLin, events] = simulateSession(scenarios(i).occMap, scenarios(i).pEventGivenX, nFrames);
    sim = analyzeSession(posLin, events, nx, ny, nShuffles);

    scenarios(i).sim = sim;
    scenarios(i).naiveSaysPlaceCell = scenarios(i).peakConditional >= naivePeakThreshold;
    scenarios(i).miSaysPlaceCell = (sum(events) >= 5) && (sim.si > sim.shuffle95);
end

%% Figure 1: occupancy, conditional probability, and local SI contribution
mapFig = figure('Color', 'w', 'Position', [80 80 1600 980]);

occStack = cat(3, scenarios.occMap);
localInfoStack = cat(3, scenarios.localInfoMap);
maxOcc = max(occStack(:));
maxLocalInfo = max(localInfoStack(:));

for i = 1:nScenarios
    subplot(3, nScenarios, i);
    imagesc(scenarios(i).occMap, [0 maxOcc]);
    axis image;
    set(gca, 'YDir', 'normal');
    title(sprintf('%s\n%s', scenarios(i).name, scenarios(i).description), ...
        'FontSize', 12);
    xlabel('Track position (x bin)');
    ylabel('Track width (y)');
    colorbar;
    colormap(gca, parula);

    subplot(3, nScenarios, i + nScenarios);
    imagesc(scenarios(i).pEventGivenX, [0 peakTarget]);
    axis image;
    set(gca, 'YDir', 'normal');
    title(sprintf('P(event|x), peak = %.2f; P(x_{peak}) = %.4f', ...
        scenarios(i).peakConditional, scenarios(i).peakBinOccupancy));
    xlabel('Track position (x bin)');
    ylabel('Track width (y)');
    colorbar;
    colormap(gca, parula);

    subplot(3, nScenarios, i + 2 * nScenarios);
    imagesc(scenarios(i).localInfoMap, [0 maxLocalInfo]);
    axis image;
    set(gca, 'YDir', 'normal');
    title(sprintf('Local SI; total = %.4f bits', scenarios(i).expectedSI));
    xlabel('Track position (x bin)');
    ylabel('Track width (y)');
    colorbar;
    colormap(gca, hot);
end

saveas(mapFig, fullfile(outDir, 'spatial_mi_linear_track_maps.png'));

%% Figure 2: along-track profiles plus summary classification
summaryFig = figure('Color', 'w', 'Position', [80 80 1550 820]);
lineColors = [0.20 0.45 0.80; 0.85 0.45 0.15; 0.70 0.20 0.20];

subplot(2, 2, 1);
hold on;
for i = 1:nScenarios
    plot(scenarios(i).trackOccupancyMarginal, 'Color', lineColors(i, :), 'LineWidth', 2);
end
hold off;
xlim([1 nx]);
xlabel('Track position (x bin)');
ylabel('Marginal occupancy');
title('Occupancy along the track');
legend({scenarios.shortName}, 'Location', 'northeast');

subplot(2, 2, 2);
hold on;
for i = 1:nScenarios
    plot(scenarios(i).trackEventMarginal, 'Color', lineColors(i, :), 'LineWidth', 2);
end
hold off;
xlim([1 nx]);
ylim([0 0.7]);
xlabel('Track position (x bin)');
ylabel('Width-averaged P(event|x)');
title('Conditional event probability along the track');

subplot(2, 2, 3);
hold on;
for i = 1:nScenarios
    plot(scenarios(i).trackLocalSIMarginal, 'Color', lineColors(i, :), 'LineWidth', 2);
end
hold off;
xlim([1 nx]);
xlabel('Track position (x bin)');
ylabel('Summed local SI contribution');
title('Where SI actually comes from');

subplot(2, 2, 4);
peakVals = [scenarios.peakConditional];
expectedVals = [scenarios.expectedSI];
simSI = arrayfun(@(s) s.sim.si, scenarios);
simThr = arrayfun(@(s) s.sim.shuffle95, scenarios);

yyaxis left;
b1 = bar(1:nScenarios, peakVals, 0.35, 'FaceColor', [0.75 0.75 0.75]);
hold on;
plot([0.5, nScenarios + 0.5], [naivePeakThreshold, naivePeakThreshold], '--k', 'LineWidth', 1.2);
ylabel('Peak P(event|x)');
ylim([0 0.7]);

yyaxis right;
b2 = bar((1:nScenarios) + 0.35, expectedVals, 0.35, 'FaceColor', [0.30 0.55 0.85]);
plot((1:nScenarios) + 0.35, simThr, 'kx', 'MarkerSize', 11, 'LineWidth', 2);
plot((1:nScenarios) + 0.35, simSI, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
ylabel('Spatial information (bits)');

for i = 1:nScenarios
    if scenarios(i).miSaysPlaceCell
        label = 'MI yes';
    else
        label = 'MI no';
    end
    text(i + 0.35, max(simSI(i), simThr(i)) * 1.08 + eps, label, ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

hold off;
set(gca, 'XTick', (1:nScenarios) + 0.18, 'XTickLabel', {scenarios.shortName});
xtickangle(20);
title('Naive peak rule versus SI');
legend([b1, b2], {'Peak P(event|x)', 'Expected SI'}, 'Location', 'northwest');

saveas(summaryFig, fullfile(outDir, 'spatial_mi_linear_track_summary.png'));

%% Command-window summary
fprintf('\nSpatial mutual information demo: nearly 2D linear track\n');
fprintf('=========================================================\n');
fprintf('Saved:\n');
fprintf('  %s\n', fullfile(outDir, 'spatial_mi_linear_track_maps.png'));
fprintf('  %s\n', fullfile(outDir, 'spatial_mi_linear_track_summary.png'));
fprintf('\n');
fprintf('%-20s  %-10s  %-10s  %-10s  %-12s  %-12s  %-12s  %-10s\n', ...
    'Scenario', 'PeakP', 'PeakOcc', 'GlobalP', 'ExpectedSI', 'EmpiricalSI', 'Shuffle95', 'MI class');
for i = 1:nScenarios
    if scenarios(i).miSaysPlaceCell
        miClass = 'yes';
    else
        miClass = 'no';
    end
    fprintf('%-20s  %-10.3f  %-10.4f  %-10.3f  %-12.4f  %-12.4f  %-12.4f  %-10s\n', ...
        scenarios(i).shortName, ...
        scenarios(i).peakConditional, ...
        scenarios(i).peakBinOccupancy, ...
        scenarios(i).globalEventProb, ...
        scenarios(i).expectedSI, ...
        scenarios(i).sim.si, ...
        scenarios(i).sim.shuffle95, ...
        miClass);
end

%% Local functions
function z = trackAlignedField(x, y, centerline, muX, sx, sy)
z = exp(-((x - muX) .^ 2) ./ (2 * sx ^ 2)) .* ...
    exp(-((y - centerline) .^ 2) ./ (2 * sy ^ 2));
end

function z = normalizeProbMap(z)
z = z ./ sum(z(:));
end

function z = maxNormalize(z)
z = z ./ max(z(:));
end

function [si, localInfoMap, pEventGlobal] = spatialMIFromMaps(pXMap, pEventGivenXMap)
pX = pXMap(:);
pEventGivenX = pEventGivenXMap(:);

pEventGlobal = sum(pX .* pEventGivenX);
pEventGlobal = min(max(pEventGlobal, eps), 1 - eps);

iPos = safeTerm(pEventGivenX, pEventGlobal) + ...
       safeTerm(1 - pEventGivenX, 1 - pEventGlobal);

localInfo = pX .* iPos;
si = sum(localInfo);
localInfoMap = reshape(localInfo, size(pXMap));
end

function t = safeTerm(q, p)
t = zeros(size(q));
mask = q > 0;
t(mask) = q(mask) .* log2(q(mask) ./ p);
end

function profile = weightedTrackAverage(valueMap, occMap)
occByX = sum(occMap, 1);
weightedSum = sum(valueMap .* occMap, 1);
profile = zeros(1, size(valueMap, 2));
valid = occByX > 0;
profile(valid) = weightedSum(valid) ./ occByX(valid);
end

function [posLin, events] = simulateSession(occMap, pEventGivenXMap, nFrames)
occVec = occMap(:);
edges = [0; cumsum(occVec)];
edges(end) = 1;

posLin = discretize(rand(nFrames, 1), edges);
pEvent = pEventGivenXMap(posLin);
events = rand(nFrames, 1) < pEvent;
end

function sim = analyzeSession(posLin, events, nx, ny, nShuffles)
nFramesLocal = numel(events);
nBins = nx * ny;

occCounts = accumarray(posLin, 1, [nBins 1], @sum, 0);
if any(events > 0)
    eventCounts = accumarray(posLin(events > 0), 1, [nBins 1], @sum, 0);
else
    eventCounts = zeros(nBins, 1);
end

pX = occCounts ./ nFramesLocal;
pEventGivenX = zeros(nBins, 1);
occupied = occCounts > 0;
pEventGivenX(occupied) = eventCounts(occupied) ./ occCounts(occupied);

[si, localInfoMap] = spatialMIFromMaps(reshape(pX, ny, nx), reshape(pEventGivenX, ny, nx));

shuffleSI = zeros(nShuffles, 1);
for s = 1:nShuffles
    shift = randi(nFramesLocal - 1, 1);
    shuffledEvents = circshift(events, shift);
    if any(shuffledEvents > 0)
        shuffledEventCounts = accumarray(posLin(shuffledEvents > 0), 1, [nBins 1], @sum, 0);
    else
        shuffledEventCounts = zeros(nBins, 1);
    end

    pEventGivenXShuffled = zeros(nBins, 1);
    pEventGivenXShuffled(occupied) = shuffledEventCounts(occupied) ./ occCounts(occupied);

    shuffleSI(s) = spatialMIFromMaps(reshape(pX, ny, nx), ...
        reshape(pEventGivenXShuffled, ny, nx));
end

shuffle95 = percentile(shuffleSI, 95);

sim = struct();
sim.si = si;
sim.shuffle95 = shuffle95;
sim.shuffleSI = shuffleSI;
sim.occCounts = reshape(occCounts, ny, nx);
sim.pEventGivenX = reshape(pEventGivenX, ny, nx);
sim.localInfoMap = localInfoMap;
sim.totalEvents = sum(events);
end

function p = percentile(values, whichPct)
sortedVals = sort(values(:));
idx = max(1, ceil((whichPct / 100) * numel(sortedVals)));
p = sortedVals(idx);
end