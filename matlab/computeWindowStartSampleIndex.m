function [windowStart, nSampWin] = computeWindowStartSampleIndex(Fs,npts,windowMin,overlapPercent)
%
% USAGE: [windowStart, nSampWin] = computeWindowStartSampleIndex(Fs,npts,windowMin,overlapPercent)
%
% Compute the start time of each Xcor-window (in samples) based on the
% correlation window length and percentage of overlap between windows
%
% INPUT:
%   Fs             = sample rate (Hz)
%   npts           = number of points total in the data
%   windowMin      = length of the correlation window in minutes
%   overlapPercent = percentage of overlap in sliding window
% OUTPUT:
%   windowStart = vector with the starting sample of each window
%   nSampWin    = number of samples in Xcor-window
%
% Written by Dylan Mikesell (mikesell@mit.edu)
% Last modified 30 June 2014

nSampWin    = windowMin * 60 * Fs; % number of sample in the window
nSlideWin   = floor(nSampWin*overlapPercent); % number of samples to move from window to next
windowStart = 1 : nSlideWin : npts ; % starting index of windows

return
