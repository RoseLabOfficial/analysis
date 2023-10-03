% Define the directory containing the files
function renaming(path, oldString, newString)
files = dir(path);

% Loop through each file
for i = 1:length(files)
    % Check if the file is not a directory and contains the substring '-0.04'
    if ~files(i).isdir && contains(files(i).name, oldString)
        % Generate the old and new filenames
        oldName = fullfile(folderPath, files(i).name);
        newName = fullfile(folderPath, strrep(files(i).name, oldString, newString));
        
        % Rename the file
        movefile(oldName, newName);
        fprintf('Renamed %s to %s\n', oldName, newName);
    end
end
fprintf('Finished renaming files.\n');
end