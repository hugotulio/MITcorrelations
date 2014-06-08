function Wout = RemoveSpikesWaveform2(W,Window,Nstd)
%
% This is a wrapper function for Piero's despiking code that works on
% single traces.
%
% USAGE: Wout = RemoveSpikesWaveform(W,catalog,Fs,Window,Nstd,FB,JulDay)
%
% INPUT:
%
% OUTPUT:
%
% Written by Dylan Mikesell (mikesell@mit.edu)
% Last modified 2 June 2014

nW   = numel(W);
Wout = W; % copy header information

for ii = 1:nW

    trace    = RemoveSpikes2(double(W(ii)), get(W(ii),'FREQ'), Window, Nstd); % despike single trace
    % add a history comment
    Wout(ii) = addhistory(Wout(ii), 'Removed spikes from data using Nstd=%0.2f and Window=%0.2f.', Nstd, Window); 
    Wout(ii) = set( Wout(ii), 'DATA', trace); % replace old data with new
end

return