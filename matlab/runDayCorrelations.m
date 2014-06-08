function  C = runDayCorrelations(W,winLength,overlapPercent)
%
% This function computes crosscorrelations of all traces in a waveform
% object at the given time window size (winLength). The windows will
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
% Last modified 2 June 2014

nW = numel(W); % number of waveforms
nC = sum(1:nW); % number of correlation pairs

C = waveform(); % blank waveform object
C = repmat(C,nC,1); % allocate complete waveform object

cnt = 0; % a counter

% double loop to cover all pairs of correlations
for ii = 1:nW
    for jj = ii:nW
        
        cnt    = cnt + 1; % update counter
        C(cnt) = correlateTwoWaveforms(W(ii),W(jj),winLength,overlapPercent);
        
    end
end

return


