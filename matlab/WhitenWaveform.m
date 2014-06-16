function Wout = WhitenWaveform(W,F)
% 
% This is a wrapper function to whiten a single waveform trace based
% on Piero's BlanchMat.m routine that worked to whiten an entire
% matrix.
%
% USAGE: Wout = WhitenWaveform(W,F)
%
% INPUT:
%   W = waveform object
%   F = frequency vector over which to whiten data [fmin,fmax] in Hz.
% Output:
%   Wout = whitened waveform object
%
% Written by Dylan Mikesell (mikesell@mit.edu)
% Last modified 2 June 2014

nW   = numel(W);
Wout = W; % allocate for speed

for ii = 1:nW
    trace    = whitenTraceWithNoZeros( double(W(ii)), F, 1/get(W(ii),'FREQ') );
    Wout(ii) = addhistory( Wout(ii), 'Whitened data over band [%0.2f, %0.2f] (Hz)', F(1), F(2) ); % add a history comment
    Wout(ii) = set( Wout(ii), 'DATA', trace); % replace old data with new
    Wout(ii) = addfield( Wout(ii), 'isWhite', true); % indicate that whitening has occurred
    Wout(ii) = addfield( Wout(ii), 'WhiteBand', F); % add whitening band parameters
end

return


