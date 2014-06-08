function trace = RemoveSpikes2(trace,Fs,Window,Nstd)
%
% Function to remove spikes from a seismic trace.
%
% USAGE: trace = RemoveSpikes(inputTrace,Fs,Window,Nstd)
%
% INPUT:
%   trace  = seismic trace
%   Fs     = sample frequency (Hz)
%   Window = length of time window (min)
%   Nstd   = number of std. deviations to set SPIKE threshold
% OUPUT:
%   trace  = de-spiked seismic trace
%
% ORIGINAL BY Piero Poli Massachussets Institute of Technology
% Last modified by Dylan Mikesell; 7 June 2014

%---------------------------------------------------------------------
% 1) load single station data

npts = numel(trace); % number of samples in the trace

%---------------------------------------------------------------------
% 2) calculate energy of seismic noise using no-null remaining window

eIdx = trace.^2 > 0; % find sample indices with non-zero energy

%---------------------------------------------------------------------
% 3) compute std. deviation of non-zero energy part of trace

Mu = std( trace(eIdx).^2 ); % get std for non-null trace (daily value of noise energy...)

spikeLimit = Nstd * Mu; % compute the signal amplitude limit that constitutes a spike

%---------------------------------------------------------------------
% 3) create windows through the data

nSampWin    = Window * 60 * Fs; % number of sample in the window
windowStart = 1 : nSampWin : npts ; % starting index of windows
% windowARray = 1 : nSampWin : (npts + 1 - nSampWin) ; % number of samples in 24 hours

nWindows = floor(npts/nSampWin); % number of windows

Tw = tukeywin( nSampWin, 0.1); % create hanning taper with 5% on front and back

%---------------------------------------------------------------------
% 7) chek if window has max energy larger than N times the daily energy STD

for ii = 1 : nWindows
    
    winSampIdx = windowStart(ii) : windowStart(ii) + nSampWin - 1; % smaple indices for this window
    tmpTrace   = trace(winSampIdx); % extract window of trace
    EnergyMax  = max(tmpTrace.^2); % get max of window
   
    if EnergyMax > spikeLimit % check if max energy greater than N*std
        
        if ii == 1 % if the first window, set window to zero
            trace(winSampIdx) = 0;
        else % else set front and back windows to zero
            trace( windowStart(ii) - nSampWin + 1 : windowStart(ii) + nSampWin - 1 ) = 0;
        end
        
    else % apply taper to edges of the window

        trace(winSampIdx) = trace(winSampIdx).* Tw; % apply taper
        
        % check how many sample in the window are larger than zero
        indx2 = abs(trace(winSampIdx)) > eps ; % D.M. replaced 1e-40 with eps (2 June 2014)
        
        % remove windows with more than 30% percent zero
        if sum(indx2) < (nSampWin - nSampWin * 0.3)
            trace(winSampIdx) = 0;
        end
        
    end
end

return