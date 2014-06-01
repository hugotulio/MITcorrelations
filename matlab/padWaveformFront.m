function Wadd = padWaveformFrontBack(Wadd)

% This function checks whether or not time series start within 1 minute of
% the start of the day and pads with zeros to one minute if not.
%
% Written by Dylan Mikesell (mikesell.at.mit.dot.edu)
% Last modified 7 May 2014
%

%%    
%--------------------------------------------------------------------------
% Check front of 

wFreq                = get(Wadd,'freq'); % sample rate
startTimeAdd         = gettimerange(Wadd); % MATLAB date integer
[~,~,~,hour,min,sec] = datevec(startTimeAdd);
secLimit             = hour*3600 + min*60 + sec; % total second starting after hour zero

if secLimit > 60
    
    D     = double(Wadd);     % get time series data
    % check row or column vector
    isflip = 0;
    if size(D,2) > 1
        D = transpose(D); % make column vector
        isflip = 1; % a flag to flip back
    end
    %--------------------------------------------------------------------------
    % build 1 minute taper make sure edges are zero
    nTap                          = floor(60*wFreq); % 60 seconds
    taper                         = tukeywin(2*nTap-1,1);
    D(1:nTap)                     = taper(1:nTap).*D(1:nTap); % taper front
    D(size(D,1)-nTap+1:size(D,1)) = taper(nTap:end).*D(size(D,1)-nTap+1:size(D,1)); % taper back
    
    
    nPad = floor(secLimit*wFreq); % number of samples to pad to get to t=0
    
end




