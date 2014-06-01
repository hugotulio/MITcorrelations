function Wadd = padWaveformFrontBack(Wadd)

% This function checks whether or not time series start within 1 minute of
% the start of the day and pads with zeros to one minute if not. It also
% checks the back of the trace to make sure it goes to the end of the day
%
% Written by Dylan Mikesell (mikesell.at.mit.dot.edu)
% Last modified 7 May 2014
%

%%    

[startTimeAdd,endTimeAdd]     = gettimerange(Wadd); % MATLAB date index
[year,month,day,hour,min,sec] = datevec(startTimeAdd);

%--------------------------------------------------------------------------
% Check front of the trace
secLimit = hour*3600 + min*60 + sec; % total second starting after hour zero
if secLimit > 60
    
    fprintf('Padding the front %0.3f seconds of waveform\n',secLimit);

    wFreq = get(Wadd,'freq'); % sample rate
    D     = double(Wadd);     % get time series data
    
    % check row or column vector
    isflip = 0;
    if size(D,2) > 1
        D      = transpose(D); % make column vector
        isflip = 1; % a flag to flip back
    end
    
    % build 1 minute taper to make sure front edge goes to zero
    nTap      = floor(60*wFreq); % 60 seconds
    taper     = tukeywin(2*nTap-1,1);
    D(1:nTap) = taper(1:nTap).*D(1:nTap); % taper front
    
    % pad the front of the trace
    nPad = floor(secLimit*wFreq); % number of samples to pad to get to t=0
    D    = padarray(D,nPad,0,'pre'); % pad front of trace
    
    % replace data and set new start time in header
    if isflip % make correct shape
        D = transpose(D);
    end
    startTime = startTimeAdd - datenum(0,0,0,0,0,secLimit);
    Wadd = set(Wadd,'data',D,'start',startTime); % set header and data
    
end

%--------------------------------------------------------------------------
% Check back of the trace

[year,month,day,hour,min,sec] = datevec(endTimeAdd);
tEnd                          = datenum(year,month,day+1,0,1,0); % should be 1 minute in the next day
[year,month,day,hour,min,sec] = datevec(tEnd-endTimeAdd);

secLimit = hour*3600 + min*60 + sec; % total second before next day
if secLimit > 60
    
    fprintf('Padding the back %0.3f seconds of waveform\n',secLimit);
    
    wFreq = get(Wadd,'freq'); % sample rate
    D     = double(Wadd);     % get time series data
    
    % check row or column vector
    isflip = 0;
    if size(D,2) > 1
        D      = transpose(D); % make column vector
        isflip = 1; % a flag to flip back
    end
    npts = size(D,1); % number of points in trace

    % build 1 minute taper to make sure front edge goes to zero
    nTap                          = floor(60*wFreq); % 60 seconds
    taper                         = tukeywin(2*nTap-1,1);
    D(size(D,1)-nTap+1:size(D,1)) = taper(nTap:end).*D(size(D,1)-nTap+1:size(D,1)); % taper back

     % replace data and set new start time in header
    if isflip % make correct shape
        D = transpose(D);
    end
    Wadd = set(Wadd,'data',D); % put tapered data back into in waveform
    % pad the end of the trace with waveform/set.m
    nPad = floor(secLimit*wFreq); % number of samples to pad to get to t=0  
    Wadd = set(Wadd,'SAMPLE_LENGTH',npts+nPad); % pad data to end

end

return

