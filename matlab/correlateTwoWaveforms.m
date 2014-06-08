% function C = correlateTwoWaveforms(W1,W2,winLength,overlapPercent)
%
% This function computes the correlation between two waveforms. The
% data are first cut based on winLength and overlapPercent. Then
% correlation occurs and the data are linearly stacked.

C = waveform(); % create blank waveform

npts1 = get(W1,'Data_Length');
npts2 = get(W2,'Data_Length');

if npts1 ~= npts2
    disp('Data lengths are the not the same. Skipping correlation.');
end

d1 = double(W1); % data vector
d2 = double(W2); % data vector

dt = 1/get(W1,'FREQ'); % sampel interval in seconds

nSampWin = floor(winLength/dt); % number of samples in a window
nOverlap = nSampWin*overlapPercent; % number of samples to move at each window

% start first window at beginning of trace
stopIdx = nSampWin;
ii = 0; % a counter

while stopIdx < npts
    
    % define sample window
    startIdx = ii*nOverlap;
    stopIdx  = ii*nOverlap + nSampWin;
    % extract samples in this window
    t1 = d1(startIdx:stopIdx); % first trace
    t2 = d2(startIdx:stopIdx); % second trace
    
    [c1,c2,c3] = normalizedCorrelation(s1, s2, Fs, SmoothMethod, Wn, K);
    % c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
    % c2: simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
    % c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
    % c4: Transfer function array normalization C12(w)/({abs(Sarray(w))^2})
    
    
    ii = ii + 1; % update counter
end

