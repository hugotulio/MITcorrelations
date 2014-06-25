function runDayCorrelations(W, windowMin, overlapPercent, smoothMethod, Wn, K, outputDirectory)
%
% This function computes crosscorrelations of all traces in a waveform
% object for the given time window size (windownMin). The windows will
% overlap based on overlapPercent value (e.g., 0.25 for 25%, etc.).
%
% USAGE: C = runDayCorrelations(W, windowMin, overlapPercent, smoothMethod, Wn, K, outputDirectory)
%
% INPUT:
%   W              = input waveform object
%   windowMin      = window length in minutes
%   overlapPercent = percent that the windows overlap
%   smoothMethod   = 'taper' or 'median'
%   Wn             = for 'taper' method, Wn=first 2*Wn discrete
%   prolate spheroidal sequences; for 'median' method, Wn=Wn order in
%   medfilt1.m
%   K              = the K most band-limited discrete prolate
%   spheroidal sequences when using the 'taper' mehtod. Default is
%   K=2*Wn-1
% OUTPUT:
%   C = a waveform objection containing all of the stacked
%   correlations with relevant META-data.
%
% Written by Dylan Mikesell (mikesell@mit.edu)
% Last modified 8 June 2014

nW = numel(W); % number of waveforms
nC = sum(1:nW); % Max number of correlation pairs per window

% set up windowing parameters
Fs          = get(W(1),'FREQ');
npts        = get(W(1),'Data_Length');
nSampWin    = windowMin * 60 * Fs; % number of sample in the window
nSlideWin   = floor(nSampWin*overlapPercent); % number of samples to move from window to next
windowStart = 1 : nSlideWin : npts ; % starting index of windows
nWindows    = numel(windowStart); % number of windows

if overlapPercent < 1 % correct index when using overlapping windows
    nWindows = nWindows - 1;
end
fprintf('Number of windows %d.\n\n',nWindows);

% loop over time windows
for tt = 1 : nWindows
    
    fprintf('Correlating window %d of %d\n',tt,nWindows);
    
    winSampIdx = windowStart(tt) : windowStart(tt) + nSampWin - 1; % smaple indices for this window
    
    % double loop to cover all pairs of correlations
    for ii = 1 %: nW
        
        WA = double(W(ii));
        % isWhitend = isfield(W(ii),'isWhite'); % check to see if data have been spectrally whitened already
        isWhitend = 0; % assume no whitening so that C2 and C3 are not computed
        
        % check that trace has less than 75% zero before doing
        % correlation
        zeroIdx = (WA(winSampIdx) == 0);
        if sum(zeroIdx) < nSampWin*0.75
            
            % loop over stations
            for jj = ii% : nW
                
                WB = double(W(jj));
                
                % check that trace has less than 75% zero before doing
                % correlation
                zeroIdx = (WB(winSampIdx) == 0);
                if sum(zeroIdx) < nSampWin*0.75
                    
                    statC = waveform(); % blank waveform object to store new correlations
                      
                    CC = normalizedCorrelation(WA(winSampIdx), WB(winSampIdx), Fs, smoothMethod, Wn, K, isWhitend);
                    % CC.c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
                    % CC.c2: simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
                    % CC.c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
                    
                    % set the basic WAVEFORM properties
                    statC = set(statC, 'FREQ', Fs);
                    statC = set(statC, 'Data_Length', numel(winSampIdx));
                    statC = set(statC, 'Station', [get(W(ii),'Station') '-' get(W(jj),'Station')]);
                    statC = set(statC, 'Channel', [get(W(ii),'Channel') '-' get(W(jj),'Channel')]);
                    statC = set(statC, 'Start', get(W(ii),'Start') + datenum(0,0,0,0,0,(windowStart(tt)-1)/Fs));
                    statC = set(statC, 'Network', [get(W(ii),'Network') '-' get(W(jj),'Network')]);
                    statC = set(statC, 'Location', [get(W(ii),'Location') '-' get(W(jj),'Location')]);
                    
                    % add station location information
                    statC = addfield(statC, 'WALA', get(W(ii),'STLA'));
                    statC = addfield(statC, 'WALO', get(W(ii),'STLO'));
                    statC = addfield(statC, 'WAEL', get(W(ii),'STEL'));
                    statC = addfield(statC, 'WBLA', get(W(jj),'STLA'));
                    statC = addfield(statC, 'WBLO', get(W(jj),'STLO'));
                    statC = addfield(statC, 'WBEL', get(W(jj),'STEL'));
                    
                    % add correlation functions and information
                    if isfield(CC,'c1')
                        statC = addfield(statC, 'c1', CC.c1);
                    end
                    if isfield(CC,'c2')
                        statC = addfield(statC, 'c2', CC.c2);
                    end
                    if isfield(CC,'c3')
                        statC = addfield(statC, 'c3', CC.c3);
                    end
                    statC = addfield(statC, 'smoothMethod', smoothMethod);
                    statC = addfield(statC, 'Wn', Wn);
                    statC = addfield(statC, 'K', K);
                    
                    station   = get(statC,'station'); % get station pair
                    startDate = [datestr(get(statC,'start'),'yyyy_mm_dd') '_window_' num2str(tt,'%03d')];
                    
                    if exist([outputDirectory '/' station],'dir') == 0
                        mkdir([outputDirectory '/' station]);
                    end              
                    
                    fname = [outputDirectory '/' station '/' startDate '.mat'];
                    save(fname,'statC','-v7.3'); % write NEW matrix file out
        
                end % end if
            end % jj loop over WB
        end % end if
    end % ii loop over WA
    
%     % write this time window output
%     fname = [outputDirectory '/' datestr(get(W(1),'start'),'YYYY_MM_DD') '_window' num2str(tt) '.mat'];
%     save(fname,'Cout');
end

return