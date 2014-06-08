function trace = RemoveSpikes(inputTrace,catalog,Fs,Window,Nstd,FB,JulDay)
%
% Function to remove EQs signal from seismic time series.
%
% USAGE: trace = RemoveSpikes(inputTrace,catalog,Fs,Window,Nstd,FB,JulDay)
%
% PROCESSING:
% 1) load single station data
% 2) filtering
% 3) windowing data
% 4) remove eqs using catalog
% 5) calculate energy of seismic noise using no-null remaining window
% 6) evaluate the STD of the no-null daily energy (should be a level of noise energy for that day)
% 7) chek if window has max energy larger than N times the daily energy STD
% 8) cosine taper the remaining no-null windows
%
%
% PARAMETERS:
% -Window: size of the window in minute.
%
% -catalog: char value of the EQ catalog used. Function to read catalog is
% at the end of this fucntion. Catlog format is: 200701010031A CARLSBERG
% RIDGE  (eq name GMT output http://www.globalcmt.org/CMTsearch.html)
%
% -Fs: sample frequency
%
% -Nstd: integer for defining std threshold
%
% CREATED BY Piero Poli Massachussets Institute of Technology
% V.1.1 - 19 Sept 2013

%---------------------------------------------------------------------
% 1) load single station data

trace        = inputTrace;
trace(end+1) = 0; % add a sample to ensure the windowing being right. Sample be removed at the end.

%---------------------------------------------------------------------
% 2) Butterworth filter the trace

[a1,b1]   = butter(2,FB*2/Fs); % filter coefficients
traceFilt = filtfilt(a1,b1,trace); % highpass filter

%---------------------------------------------------------------------
% 3) windowing data

% 3.1 define the windows
Win   = Window * 60 * Fs;
t0win = 1 : Win : (86400*Fs+1 - Win) ;

% 3.2 Get information from seismic catalog
out = readCMTCatalog(catalog,'no',2007);

% 3.2.1 find earthquakes in this day
A = find( out.julianD == JulDay );

% 3.2.2 get time of earthquakes
houreqs   = (str2num(out.hour(A,:))) * 3600 * Fs; % hour of earthquake (in seconds)
minuteeqs = (str2num(out.minute(A,:))) * 60 * Fs; % minute of earthquake (in seconds)
EqTime    = houreqs + minuteeqs; % time of earthquake for this day (in seconds)

clear A

%---------------------------------------------------------------------
% 4) window analysis of eqs presence
for i2 = 1 : numel(t0win) - 1
    
    indx = find( EqTime > t0win(i2) && EqTime < t0win(i2+1), 1);
    
    % remove eqs using catalog
    if isempty(indx) == 0
        if i2 > 21
            trace(t0win(i2):end) = 0;
        else
            trace(t0win(i2):t0win(i2+2)) = 0;
        end
    end
end

%---------------------------------------------------------------------
% 5) calculate energy of seismic noise using no-null remaining window

ee = traceFilt.^2 > 0; % find samples with non-zero energy

%---------------------------------------------------------------------
% 6) evaluate the STD of the no-null daily energy (should be a level of noise energy for that day)

Mu = std( traceFilt(ee).^2 ); % get std for non-null trace (daily value of noise energy...)

%---------------------------------------------------------------------
% 7) chek if window has max energy larger than N times the daily energy STD

Tw   = tukeywin(Win+1,.1); % allocate hanning taper
wabs = zeros( numel(t0win), 1 ); % allocate for speed

for i3 = 1 : numel(t0win) - 1
    
    win = traceFilt( t0win(i3) : (t0win(i3)+Win) );
    
    wabs(i3) = max( (win).^2 ); % get max of window
    
    winSampIdx = t0win(i3) : t0win(i3) + Win;
    
    if wabs(i3) > Nstd * Mu % check if max greater than N*std
        
        if i3 == 1
            trace(winSampIdx) = 0;
        else
            trace( t0win(i3) - Win : t0win(i3) + Win ) = 0;
        end
        
    else
  
        trace(winSampIdx) = trace(winSampIdx).* Tw; % apply taper
        
        % check how many sample in the window are larger than zero
        indx2 = find( abs(trace(winSampIdx)) > eps ); % D.M. replaced 1e-40 with eps (2 June 2014)
        
        % remove windows with zeros
        if numel(indx2) < Win - 1000
            trace(winSampIdx) = 0;
        end
        
    end
    clear win
end

trace(end)=[]; % remove the sample added at the beginning...

return