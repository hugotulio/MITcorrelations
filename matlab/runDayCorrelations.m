function  C = runDayCorrelations(W,windowMin,overlapPercent,smoothMethod,Wn,K,outputDirectory)
%
% This function computes crosscorrelations of all traces in a waveform
% object for the given time window size (winLength). The windows will
% overlap based on overlapPercent. The correlations are linearly
% stacked over the day-long waveform traces. More complex stacking
% strategies can be applied to the day-long correlations later.
%
% USAGE: Wout = runCorrelations(W,winLength,overlapPercent)
%
% INPUT:
%   W              = input waveform object, assuming 1 day length
%   winLength      = window length in seconds
%   overlapPercent = percent that the windows overlap
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
nSlideWin   = floor(nSampWin*overlapPercent); % number of samples to move from window to nex
windowStart = 1 : nSlideWin : npts ; % starting index of windows
nWindows    = numel(windowStart) - 1; % number of windows

% loop over time windows
for tt = 1:nWindows
    
    winSampIdx = windowStart(tt) : windowStart(tt) + nSampWin - 1; % smaple indices for this window
    
    cnt = 0; % a counter
    C = waveform(); % blank waveform object to store new correlations
    C = repmat(C,nC,1); % allocate complete waveform object

    % double loop to cover all pairs of correlations
    for ii = 1:nW
        
        WA = double(W(ii));
        
        for jj = ii:nW
            cnt    = cnt + 1; % update counter
            
            WB = double(W(jj));
            
            [c1,c2,c3] = normalizedCorrelation(WA(winSampIdx), WB(winSampIdx), Fs, smoothMethod, Wn, K);
            % c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
            % c2: simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
            % c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
            % c4: Transfer function array normalization C12(w)/({abs(Sarray(w))^2})
            
            % set the basic WAVEFORM properties
            C(cnt) = set(C(cnt), 'FREQ', Fs);
            C(cnt) = set(C(cnt), 'Data_Length', numel(winSampIdx));
            C(cnt) = set(C(cnt), 'Station', [get(W(ii),'Station') '-' get(W(jj),'Station')]);
            C(cnt) = set(C(cnt), 'Channel', [get(W(ii),'Channel') '-' get(W(jj),'Channel')]);
            C(cnt) = set(C(cnt), 'Start', get(W(ii),'Start') + datenum(0,0,0,0,0,(windowStart(tt)-1)/Fs));
            C(cnt) = set(C(cnt), 'Network', [get(W(ii),'Network') '-' get(W(jj),'Network')]);
            C(cnt) = set(C(cnt), 'Location', [get(W(ii),'Location') '-' get(W(jj),'Location')]);
  
            % add station location information
            C(cnt) = addfield(C(cnt), 'WALA', get(W(ii),'STLA'));
            C(cnt) = addfield(C(cnt), 'WALO', get(W(ii),'STLO'));
            C(cnt) = addfield(C(cnt), 'WAEL', get(W(ii),'STEL'));
            C(cnt) = addfield(C(cnt), 'WBLA', get(W(jj),'STLA'));
            C(cnt) = addfield(C(cnt), 'WBLO', get(W(jj),'STLO'));
            C(cnt) = addfield(C(cnt), 'WBEL', get(W(jj),'STEL'));

            % add correlation information
            C(cnt) = addfield(C(cnt), 'c1', c1);
            C(cnt) = addfield(C(cnt), 'c2', c2);
            C(cnt) = addfield(C(cnt), 'c3', c3);
            C(cnt) = addfield(C(cnt), 'smoothMethod', smoothMethod);
            C(cnt) = addfield(C(cnt), 'Wn', Wn);
            C(cnt) = addfield(C(cnt), 'K', K);
                       
        end
    end
    
    % write this time window output
    fname = [outputDirectory '/' datestr(get(W(1),'start'),'YYYY_MM_DD') '_window' num2str(tt) '.mat'];
    save(fname,'C');
end

return


