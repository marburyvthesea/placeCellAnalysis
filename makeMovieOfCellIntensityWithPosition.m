function makeMovieOfCellIntensityWithPosition(imageData, scalingFactors, x_position, y_position, outputFileName)
    % Number of frames
    numFrames = numel(scalingFactors);

    % Create a movie object
    movieObj = VideoWriter(outputFileName);
    open(movieObj);

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

    % Iterate through each scaling factor and create frames
    for frameIndex = 1:numFrames
        % Multiply the image matrix by the scaling factor
        scaledImage = imageData * scalingFactors(frameIndex);

        % Display the movie frame
        set(hMovie, 'CData', scaledImage);

        % Update the x and y position plot with only one dot at a time
        set(hPosition, 'XData', x_position(frameIndex), 'YData', y_position(frameIndex));

        % Capture the frame for the movie
        frame = getframe(gcf);

        % Write the frame to the movie file
        writeVideo(movieObj, frame);
    end

    % Close the movie file
    close(movieObj);
end