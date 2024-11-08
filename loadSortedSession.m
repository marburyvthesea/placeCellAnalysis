
% Open a file dialog for the user to select a CSV file
[fileName, filePath] = uigetfile('*.csv', 'Select CSV File');

% Check if the user clicked Cancel
if isequal(fileName, 0)
    disp('User canceled the operation');
else
    % Use readtable to load the selected CSV file into a table
    csvFilePath = fullfile(filePath, fileName);
    dataTable = readtable(csvFilePath, 'VariableNamesLine', 1);
    dataTable(1, :) = [];
end