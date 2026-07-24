% spatial_mutual_information_demo.m
%
% Synthetic demonstration of why occupancy-weighted spatial mutual
% information is safer than using only a high conditional firing
% probability P(event | bin) to classify place cells.
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
% The exact equations were reproduced in the Methods of:
%   Kinsky et al., Nature Communications (2020)
%   "Trajectory-modulated hippocampal neurons persist throughout
%    memory-guided navigation", doi: 10.1038/s41467-020-16226-4
%
% What this script shows:
%   1) A true place cell with a compact field and ordinary occupancy.
%   2) A fast pass-through hotspot: high speed causes very low occupancy in
%      the hotspot bin(s), so peak P(event | bin) looks impressive but SI
%      stays small because P(x) is tiny.
%   3) A slow linger zone: low speed causes high occupancy in a reward-like
%      zone, but the cell only shows weak contrast relative to its already
%      high global event rate, so SI stays small even though peak
%      P(event | bin) is high.
%
% Output:
%   - spatial_mi_demo_maps.png
%   - spatial_mi_demo_summary.png
%
% Notes:
%   - This is a pedagogical simulation, not a re-analysis of the paper.
%   - The place-cell significance test uses shuffled event timestamps,
%     mirroring the paper's idea. Here we use circular shifts.

clear;
clc;
close all;

rng(7);

%% User-tunable parameters
nx = 20;                    % spatial bins in x
ny = 20;                    % spatial bins in y
nFrames = 3000;             % synthetic imaging frames
nShuffles = 1000;           % as in the paper's shuffle logic
naivePeakThreshold = 0.50;  % intentionally naive "high P(event|bin)" rule
peakTarget = 0.60;          % all synthetic cells share this peak probability
outDir = pwd;

%% Base occupancy map P(x)
[xGrid, yGrid] = meshgrid(1:nx, 1:ny);
centerX = (nx + 1) / 2;
centerY = (ny + 1) / 2;

baseOccMap = 0.20 + 0.80 * gaussian2d(xGrid, yGrid, centerX, centerY, 4.2, 4.2);
baseOccMap = baseOccMap ./ sum(baseOccMap(:));

% Increase 0.96 to make the fast-transit hotspot even less occupied.
fastTransitOccMap = baseOccMap .* ...
    (1 - 0.96 * normalized(gaussian2d(xGrid, yGrid, 2.0, 18.5, 0.9, 0.9))) + 1e-5;
fastTransitOccMap = fastTransitOccMap ./ sum(fastTransitOccMap(:));

% Increase 1.7 to make the low-velocity linger zone even more occupied.
slowLingerOccMap = baseOccMap + 1.7 * gaussian2d(xGrid, yGrid, 15.0, 5.0, 3.0, 3.0);
slowLingerOccMap = slowLingerOccMap ./ sum(slowLingerOccMap(:));

%% Synthetic cells
placeMap = 0.02 + (peakTarget - 0.02) * ...
    normalized(gaussian2d(xGrid, yGrid, 6.5, 14.5, 2.2, 2.0));

fastTransitHotspotMap = 0.02 + (peakTarget - 0.02) * ...
    normalized(gaussian2d(xGrid, yGrid, 2.0, 18.5, 0.45, 0.45));

slowLingerStateMap = 0.45 + (peakTarget - 0.45) * ...
    normalized(gaussian2d(xGrid, yGrid, 15.0, 5.0, 3.0, 3.0));

scenarios = struct([]);
scenarios(1).name = 'True place field';
scenarios(1).shortName = 'True field';
scenarios(1).occMap = baseOccMap;
scenarios(1).pEventGivenX = placeMap;
scenarios(1).description = 'Compact field, low baseline';

scenarios(2).name = 'Fast pass-through hotspot';
scenarios(2).shortName = 'Fast transit';
scenarios(2).occMap = fastTransitOccMap;
scenarios(2).pEventGivenX = fastTransitHotspotMap;
scenarios(2).description = 'High speed produces very low occupancy in the hotspot';

scenarios(3).name = 'Slow linger zone';
scenarios(3).shortName = 'Slow linger';
scenarios(3).occMap = slowLingerOccMap;
scenarios(3).pEventGivenX = slowLingerStateMap;
scenarios(3).description = 'Low speed raises occupancy, but spatial contrast is weak';

nScenarios = numel(scenarios);

%% Expected SI from the underlying maps
for i = 1:nScenarios
    [expectedSI, localInfoMap, pEventGlobal] = spatialMIFromMaps( ...
        scenarios(i).occMap, scenarios(i).pEventGivenX);

    scenarios(i).expectedSI = expectedSI;
    scenarios(i).localInfoMap = localInfoMap;
    scenarios(i).globalEventProb = pEventGlobal;
    scenarios(i).peakConditional = max(scenarios(i).pEventGivenX(:));
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
mapFig = figure('Color', 'w', 'Position', [80 80 1500 980]);

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
    xlabel('x bin');
    ylabel('Occupancy P(x)');
    colorbar;
    colormap(gca, parula);

    subplot(3, nScenarios, i + nScenarios);
    imagesc(scenarios(i).pEventGivenX, [0 peakTarget]);
    axis image;
    set(gca, 'YDir', 'normal');
    title(sprintf('P(event|x), peak = %.2f', scenarios(i).peakConditional));
    xlabel('x bin');
    ylabel('Conditional prob.');
    colorbar;
    colormap(gca, parula);

    subplot(3, nScenarios, i + 2 * nScenarios);
    imagesc(scenarios(i).localInfoMap, [0 maxLocalInfo]);
    axis image;
    set(gca, 'YDir', 'normal');
    title(sprintf('Local SI; total = %.4f bits', scenarios(i).expectedSI));
    xlabel('x bin');
    ylabel('SI contribution');
    colorbar;
    colormap(gca, hot);
end

saveas(mapFig, fullfile(outDir, 'spatial_mi_demo_maps.png'));

%% Figure 2: summary comparison
summaryFig = figure('Color', 'w', 'Position', [80 80 1450 470]);
barColors = [0.20 0.45 0.80; 0.85 0.45 0.15; 0.70 0.20 0.20];

subplot(1, 3, 1);
peakVals = [scenarios.peakConditional];
b1 = bar(1:nScenarios, peakVals, 0.65, 'FaceColor', 'flat');
for i = 1:nScenarios
    b1.CData(i, :) = barColors(i, :);
end
hold on;
plot([0.5, nScenarios + 0.5], [naivePeakThreshold, naivePeakThreshold], '--k', ...
    'LineWidth', 1.2);
text(0.65, naivePeakThreshold + 0.015, ...
    sprintf('naive threshold = %.2f', naivePeakThreshold), ...
    'FontSize', 10, 'HorizontalAlignment', 'left');
hold off;
set(gca, 'XTick', 1:nScenarios, 'XTickLabel', {scenarios.shortName});
ylim([0 0.7]);
ylabel('Peak P(event|x)');
title('Simple conditional probability');
xtickangle(20);

subplot(1, 3, 2);
expectedVals = [scenarios.expectedSI];
b2 = bar(1:nScenarios, expectedVals, 0.65, 'FaceColor', 'flat');
for i = 1:nScenarios
    b2.CData(i, :) = barColors(i, :);
end
set(gca, 'XTick', 1:nScenarios, 'XTickLabel', {scenarios.shortName});
ylabel('Expected SI (bits)');
title('Occupancy-weighted spatial MI');
xtickangle(20);

subplot(1, 3, 3);
simSI = arrayfun(@(s) s.sim.si, scenarios);
simThr = arrayfun(@(s) s.sim.shuffle95, scenarios);
b3 = bar(1:nScenarios, simSI, 0.65, 'FaceColor', 'flat');
for i = 1:nScenarios
    b3.CData(i, :) = barColors(i, :);
end
hold on;
hThr = plot(1:nScenarios, simThr, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
for i = 1:nScenarios
    if scenarios(i).miSaysPlaceCell
        label = 'MI says yes';
    else
        label = 'MI says no';
    end
    text(i, max(simSI(i), simThr(i)) * 1.05 + eps, label, ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end
hold off;
set(gca, 'XTick', 1:nScenarios, 'XTickLabel', {scenarios.shortName});
ylabel('Empirical SI (bits)');
title('One simulated session: SI vs. 95th shuffle percentile');
legend([b3, hThr], {'Empirical SI', 'Shuffle 95%'}, 'Location', 'northwest');
xtickangle(20);

saveas(summaryFig, fullfile(outDir, 'spatial_mi_demo_summary.png'));

%% Command-window summary
fprintf('\nSpatial mutual information demo\n');
fprintf('===============================================\n');
fprintf('Saved:\n');
fprintf('  %s\n', fullfile(outDir, 'spatial_mi_demo_maps.png'));
fprintf('  %s\n', fullfile(outDir, 'spatial_mi_demo_summary.png'));
fprintf('\n');
fprintf('%-20s  %-10s  %-10s  %-12s  %-12s  %-12s  %-10s\n', ...
    'Scenario', 'PeakP', 'GlobalP', 'ExpectedSI', 'EmpiricalSI', 'Shuffle95', 'MI class');
for i = 1:nScenarios
    if scenarios(i).miSaysPlaceCell
        miClass = 'yes';
    else
        miClass = 'no';
    end
    fprintf('%-20s  %-10.3f  %-10.3f  %-12.4f  %-12.4f  %-12.4f  %-10s\n', ...
        scenarios(i).shortName, ...
        scenarios(i).peakConditional, ...
        scenarios(i).globalEventProb, ...
        scenarios(i).expectedSI, ...
        scenarios(i).sim.si, ...
        scenarios(i).sim.shuffle95, ...
        miClass);
end

%% Local functions
function z = gaussian2d(x, y, mx, my, sx, sy)
z = exp(-(((x - mx) .^ 2) ./ (2 * sx ^ 2) + ((y - my) .^ 2) ./ (2 * sy ^ 2)));
end

function z = normalized(z)
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
