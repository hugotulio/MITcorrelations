function sortDayCorrelations3( inputDirectory, outputDirectory)
%
% This function takes time window correlation matrices in a given
% directory and sorts the time windows into new matrices based on
% station correlation pairs. In essence, we rewrite the correlation
% matrices so that we have all time windows for a given station pair
% in one file. This allows us to try different stacking later.

%---------------------------------------------------------------------
% Check input directory
if exist(inputDirectory,'dir')
    fprintf('Searching %s...\n',inputDirectory);
else
    fprintf('%s does not exist. Check PATH.\n',inputDirectory);
    return
end

%---------------------------------------------------------------------
% Check output directory
[success,message,messageID] = checkOutputDir(outputDirectory);

%---------------------------------------------------------------------
% Process matrix files
corrMatrices = dir([inputDirectory '/*.mat']); % get list of *.mat files
nMat = numel(corrMatrices);
fprintf('Found %d correlation matrices.\n',nMat);

% loop over all matrix files found
for ii = 1 : nMat
    
    % open a correlation matrix file
    fprintf('Sorting %s\n',corrMatrices(ii).name);
    load([inputDirectory '/' corrMatrices(ii).name]);
    
    % remove any empty parts first
    statPair = get(Cout,'Station');
    emptyCells = cellfun(@isempty,statPair); % find empty cells
    Cout(emptyCells) = []; % remove empty cells
    statPair(emptyCells) = [];
    
    nPairs = numel(Cout); % number of non-empty station pairs
    
    for jj = 1 : nPairs
        
        station = get(Cout(jj),'station'); % get station pair
        
        fname = [outputDirectory '/' station '.mat'];
        
        S = struct(Cout(jj));
        
        if exist(fname,'file')         
            matObj = matfile(fname, 'Writable', true); % open matfile and make writable
            nWave = numel(matObj.S); % number of waveforms in matfile already
            matObj.S(1,nWave + 1) = S; % append another waveform
        else
            save(fname, 'S', '-v7.3'); % save matfile
        end
        
    end % jj over station pairs
    
end % ii over time windows

return