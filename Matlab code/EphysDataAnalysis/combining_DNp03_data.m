% File list and corresponding variable names to extract
%For DNp03 data
fileVarMap = {
    '1st DNp03 current injection_05-21-2025.mat', {'DNp03_1of5_all_traces'};
    '2nd_3rd DNp03 current injection_05-29-2025.mat', {'DNp03_2of5_all_traces', 'DNp03_3of5_all_traces'};
    '4tth DNp03 current injection_05-30-2025.mat', {'DNp03_4of5_all_traces'};
};

% Preallocate combined data structure
combinedData = struct();

% Loop over the files
for i = 1:size(fileVarMap, 1)
    fileName = fileVarMap{i, 1};
    varList = fileVarMap{i, 2};

    % Load the file
    fileData = load(fileName);
    
    % Loop through each variable to extract
    for j = 1:length(varList)
        originalVar = varList{j};

        % Extract fly number from the variable name (e.g., 'DNp03_2of5_all_traces')
        tokens = regexp(originalVar, 'DNp03_(\d+)of5_all_traces', 'tokens');
        if isempty(tokens)
            error('Could not extract fly number from variable name: %s', originalVar);
        end
        flyNum = str2double(tokens{1}{1});

        % Rename to DNp03_flyN_all_traces
        newVarName = sprintf('Fly%d_all_traces', flyNum);

        % Add to combined structure
        if isfield(fileData, originalVar)
            combinedData.(newVarName) = fileData.(originalVar);
        else
            error('Variable %s not found in file %s.', originalVar, fileName);
        end
    end
end

% Save the combined structure
save('combined_fly_traces_DNp03.mat', '-struct', 'combinedData');
