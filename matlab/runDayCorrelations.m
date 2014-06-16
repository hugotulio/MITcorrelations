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
    
    cnt  = 0; % a counter
    Cout = waveform(); % blank waveform object to store new correlations
    Cout = repmat(Cout,nC,1); % allocate complete waveform object
    
    % double loop to cover all pairs of correlations
    for ii = 1 %: nW
        
        WA = double(W(ii));
        isWhitend = isfield(W(ii),'isWhite'); % check to see if data have been spectrally whitened already
        
        % check that trace has less than 75% zero before doing
        % correlation
        zeroIdx = (WA(winSampIdx) == 0);
        if sum(zeroIdx) < nSampWin*0.75
            
            % loop over stations
            for jj = ii : nW
                
                
                WB = double(W(jj));
                
                % check that trace has less than 75% zero before doing
                % correlation
                zeroIdx = (WB(winSampIdx) == 0);
                if sum(zeroIdx) < nSampWin*0.75
                    
                    cnt = cnt + 1; % update counter
                    
                    CC = normalizedCorrelation(WA(winSampIdx), WB(winSampIdx), Fs, smoothMethod, Wn, K, isWhitend);
                    % CC.c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
                    % CC.c2: simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
                    % CC.c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
                    
                    % set the basic WAVEFORM properties
                    Cout(cnt) = set(Cout(cnt), 'FREQ', Fs);
                    Cout(cnt) = set(Cout(cnt), 'Data_Length', numel(winSampIdx));
                    Cout(cnt) = set(Cout(cnt), 'Station', [get(W(ii),'Station') '-' get(W(jj),'Station')]);
                    Cout(cnt) = set(Cout(cnt), 'Channel', [get(W(ii),'Channel') '-' get(W(jj),'Channel')]);
                    Cout(cnt) = set(Cout(cnt), 'Start', get(W(ii),'Start') + datenum(0,0,0,0,0,(windowStart(tt)-1)/Fs));
                    Cout(cnt) = set(Cout(cnt), 'Network', [get(W(ii),'Network') '-' get(W(jj),'Network')]);
                    Cout(cnt) = set(Cout(cnt), 'Location', [get(W(ii),'Location') '-' get(W(jj),'Location')]);
                    
                    % add station location information
                    Cout(cnt) = addfield(Cout(cnt), 'WALA', get(W(ii),'STLA'));
                    Cout(cnt) = addfield(Cout(cnt), 'WALO', get(W(ii),'STLO'));
                    Cout(cnt) = addfield(Cout(cnt), 'WAEL', get(W(ii),'STEL'));
                    Cout(cnt) = addfield(Cout(cnt), 'WBLA', get(W(jj),'STLA'));
                    Cout(cnt) = addfield(Cout(cnt), 'WBLO', get(W(jj),'STLO'));
                    Cout(cnt) = addfield(Cout(cnt), 'WBEL', get(W(jj),'STEL'));
                    
                    % add correlation functions and information
                    if isfield(CC,'c1')
                        Cout(cnt) = addfield(Cout(cnt), 'c1', CC.c1);
                    end
                    if isfield(CC,'c2')
                        Cout(cnt) = addfield(Cout(cnt), 'c2', CC.c2);
                    end
                    if isfield(CC,'c3')
                        Cout(cnt) = addfield(Cout(cnt), 'c3', CC.c3);
                    end
                    Cout(cnt) = addfield(Cout(cnt), 'smoothMethod', smoothMethod);
                    Cout(cnt) = addfield(Cout(cnt), 'Wn', Wn);
                    Cout(cnt) = addfield(Cout(cnt), 'K', K);
                end % end if
            end % jj loop over WB
        end % end if
    end % ii loop over WA
    
    % write this time window output
    fname = [outputDirectory '/' datestr(get(W(1),'start'),'YYYY_MM_DD') '_window' num2str(tt) '.mat'];
    save(fname,'Cout');
end

return