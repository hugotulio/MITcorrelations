function W = addJulDayWaveform(Wadd,julDayDir,julDay,resampleFrequency)

% This function adds the waveform data in W to the Julian Day data
% Structure. Wadd is a waveform toolbox object

%**************************************************************************
% still need to correct tPad if the data start more than 1 minute after
% the hour. This might happen and now I assume all data start close to t=0
%**************************************************************************
%
% Written by Dylan Mikesell (mikesell@mit.edu)
% Last modified 7 May 2014
%
% Requires MATLAB waveform toolbox or GISMOTools
% (http://www.giseis.alaska.edu/Seis/EQ/tools/GISMO/)

%--------------------------------------------------------------------------
% Process the trace and extract
Wadd  = demean(Wadd);     % remove mean of this waveform
Wadd  = detrend(Wadd);    % remove any linear trend
Wadd  = fillgaps(Wadd,0); % fill data gaps with zeros
Wadd  = padWaveformFrontBack(Wadd); % make sure edges go beginning and end of day
%--------------------------------------------------------------------------
% check row or column vector
wFreq  = get(Wadd,'freq'); % input waveform frequency
D      = double(Wadd);     % get time series data
isflip = 0;
if size(D,2) > 1
    D      = transpose(D); % make column vector
    isflip = 1; % a flag to flip back
end
%--------------------------------------------------------------------------
% build 1 minute taper make sure edges are zero
nTap                          = floor(60*wFreq); % 60 seconds
taper                         = tukeywin(2*nTap-1,1);
D(1:nTap)                     = taper(1:nTap).*D(1:nTap); % taper front
D(size(D,1)-nTap+1:size(D,1)) = taper(nTap:end).*D(size(D,1)-nTap+1:size(D,1)); % taper back
%--------------------------------------------------------------------------
% Pad the front and back of the trace
tPad = 2*60; % (s) number of seconds to pad data
tNum = datenum(0,0,0,0,0,tPad); % matlab date number for padding length
nPad = floor(tPad*wFreq);     % number of samples to pad front and back
D    = padarray(D,nPad,0,'both'); % pad front and back of trace
%--------------------------------------------------------------------------
% set the correct times in Waveform object
startTimeAdd = gettimerange(Wadd); % MATLAB date integer
startTime    = startTimeAdd - tNum; % adjust start time
Wadd         = set(Wadd,'start',startTime); % set new start time
%--------------------------------------------------------------------------
% Resample waveform
Q    = round(wFreq/resampleFrequency); % resample factor
D    = resample(D,1,Q); % resample the data
if isflip % make correct shape
    D = transpose(D);
end
Wadd = set(Wadd,'data',D,'Freq',resampleFrequency); % set header and data
%--------------------------------------------------------------------------
% Align the data on t=00.000 for this day
[year,month,day] = datevec(startTimeAdd); % get the date of this data
startNum         = datenum(year,month,day,0,0,0); % matlab date number for aligning trace
Wadd             = align(Wadd,datestr(startNum,'mm/dd/yyyy HH:MM:SS.FFF'),get(Wadd,'freq'));
%--------------------------------------------------------------------------
% Extract the data we want for this day
startNum = datenum(year,month,day,0,0,0); % matlab date number for aligning trace
endNum   = datenum(year,month,day+1,0,0,0); % matlab date number for aligning trace
Wadd     = extract(Wadd,'time',startNum,endNum);
%--------------------------------------------------------------------------
% Add this data to the julian day matrix
fname = [julDayDir '/julDay_' num2str(julDay) '.mat'];

if exist(fname,'file')
    load(fname); % append new data in W
    wIdx = strcmp(get(W,'station'),get(Wadd,'station'));
    if sum(wIdx) == 0
        fprintf('Adding new waveform from station %s\n',get(Wadd,'station'));
        W(numel(W)+1) = Wadd;
    else
        fprintf('Station %s already in matrix\n',get(Wadd,'station'));
        fprintf('Replacing waveform with current waveform\n');
        W(wIdx) = Wadd;
    end
else
    fprintf('Adding new waveform from station %s\n',get(Wadd,'station'));
    W = Wadd; % make new waveform for output
end
%--------------------------------------------------------------------------
% write matrix
save(fname,'-v7.3','W');

return