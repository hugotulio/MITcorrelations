% function sortDayCorrelations( inputDirectory,stationListFile)
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
    corrMatrices = dir([inputDirectory '/*.mat']); % get list of *.mat files
else
    fprintf('%s does not exist. Check PATH.\n',inputDirectory);
    return
end
%---------------------------------------------------------------------
% Check station file
if exist(stationListFile,'file')
    fprintf('Loading %s...\n',stationListFile);
else
    fprintf('%s does not exist. Check PATH/FILENAME.\n',stationListFile);
    return
end   
%---------------------------------------------------------------------
% Process matrix files

nMat = numel(corrMatrices);
fprintf('Found %d correlation matrices.\n',nMat);

for ii = 1% : nMat
    load([inputDirectory '/' corrMatrices(ii).name]);
    
    nPairs = numel(C);
    
    for jj = 1 %: nPairs
        stations = get(C(jj),'Station');
        stat1 =  
    end
    
    
    
    
end