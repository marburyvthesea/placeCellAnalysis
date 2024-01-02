% load position and spike information

% m233_2023_09_29_PosSpikesAligned.mat

% load cell location information 

% 17-Nov_17_23_38_sorted_output.mat

cell_idx=48; 
imageData=A(:,:,cell_idx);

scalingFactors=C(cell_idx,:);
x_position=pos(2,:);
y_position=pos(3,:);


% run make movie function
makeMovieOfCellIntensityWithPosition(imageData, scalingFactors, x_position, y_position, 'cell_idx_movie.avi')


