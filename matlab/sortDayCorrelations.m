function sortDayCorrelations( inputDirectory, outputDirectory)
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
        
        station = get(Cout(jj),'station');
        fname = [outputDirectory '/' station '.mat'];
        
        if exist(fname,'file')
            load(fname); % load existing file
            nWin = numel(statC);
            statC( nWin + 1 ) = Cout(jj); % append this time window
        else
            statC = Cout(jj);
        end
        save(fname,'statC'); % write matrix file out
    
    end % jj over station pairs
    
end % ii over time windows

return