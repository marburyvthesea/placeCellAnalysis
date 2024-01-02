% Assuming 'imageData' is a 600x600 double
% Assuming 'scalingFactors' is a 1xN double
% Assuming 'x_position' and 'y_position' are both 1x11858 double

% Number of frames
numFrames = numel(scalingFactors);

% Create a figure for the movie and position plot
figure;

% Create a subplot for the movie
subplot(2, 1, 1);
hMovie = imshow(imageData * scalingFactors(1));  % Initialize the movie frame
title('Intensity Movie');

% Create a subplot for the x and y position plot
subplot(2, 1, 2);
hPosition = plot(x_position(1), y_position(1), 'ro');  % Initialize the position dot
axis tight;
title('Object Position');
xlabel('X Coordinate');
ylabel('Y Coordinate');

% Set constant x and y limits for the position plot
xlim([min(x_position), max(x_position)]);
ylim([min(y_position), max(y_position)]);

% Create a uifigure
uif = uifigure('Name', 'Movie and Position Plot Viewer', 'Position', [100, 100, 800, 600]);

% Create a slider for frame selection
slider = uislider(uif, 'Position', [100, 20, 600, 3], 'Limits', [1, numFrames], 'ValueChangedFcn', @(sld, event) updateFrames(sld, hMovie, hPosition, scalingFactors, x_position, y_position));

% Initialize the slider value
slider.Value = 1;

% Display the initial frame
updateFrames(slider, hMovie, hPosition, scalingFactors, x_position, y_position, imageData);

% Callback function for slider value change
function updateFrames(slider, hMovie, hPosition, scalingFactors, x_position, y_position, imageData)
    % Get the selected frame index from the slider
    frameIndex = round(slider.Value);
    
    % Update the movie frame
    set(hMovie, 'CData', imageData * scalingFactors(frameIndex));
    
    % Update the position plot
    set(hPosition, 'XData', x_position(frameIndex), 'YData', y_position(frameIndex));
end