function Wout = RemoveSpikesWaveform(W,catalog,Window,Nstd,FB,JulDay)
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
    
    Fs = get(W(ii),'FREQ'); % sample frequency
%     JulDay = we should compute from Wavefrom start time
    
    trace    = RemoveSpikes2(double(W(ii)), catalog, Fs, Window, Nstd, FB, JulDay); % despike single trace
    Wout(ii) = addhistory(Wout(ii), 'Removed spikes from data using %s catalog.', catalog); % add a history comment
    Wout(ii) = set( Wout(ii), 'DATA', trace); % replace old data with new
end

return