clear all
close all
clc

% Example to load and process the Earthquake EVENT-CATALOG

catalogName = 'data/EVENT-CATALOG';

% load the catalog
Catalog = ParserCatalog2(catalogName);

% to access the fields you use something like this

% idx = 1; % get the first event
% Catalog.number(idx)      % e.g. 0
% Catalog.catalogType(idx) % e.g. 'US'
% Catalog.magnitude(idx)   % e.g. 4
% Catalog.depth(idx)       % e.g. -11.5000
% Catalog.lat(idx)         % e.g. 38.4065
% Catalog.lon(idx)         % e.g. -119.3642
% Catalog.time1(idx);      % e.g. '20140530_0'
% Catalog.time2(idx)       % e.g. '2014-05-30T07'

%% convert date to internal matlab date format

date = char([Catalog.time2]); % convert cell to character string

year   = str2num(date(:,1:4)); % extract year from string and convert to number
mon    = str2num(date(:,6:7));
day    = str2num(date(:,9:10));
hour   = str2num(date(:,12:13));
minute = zeros(size(year,1),1); % need to make minutes and seconds for datenum.m input
second = minute;

matlabDate = datenum(year,mon,day,hour,minute,second); % create the matlab number

h = figure;
plot(matlabDate-min(matlabDate),[Catalog.magnitude],'*');
xlabel('Day (0 = 1 Jan. 2006)'); ylabel('Magnitude');

%% load a day matrix

load('data/julDay_1952.mat');

%% plot

C = correlation(W);
plot(C);

%% look at window start times

windowMin      = 60; % length in minutes of correlation window 
overlapPercent = 0.5; % percentage of overlap for windows

% set up windowing parameters and compute start time of each window
% (in samples)
[windowStart, nSampWin] = computeWindowStartSampleIndex(get(W(1),'FREQ'),get(W(1),'Data_Length'),windowMin,overlapPercent);

% absolute start times for all correlation windows
absStartTimes = get(W(1),'Start') + datenum(0,0,0,0,0,(windowStart/Fs));

%% compare matlabDate and absStartTimes

% This way will figure out which correlation windows contain earthquakes
clc

windowMatlab = datenum(0,0,0,0,windowMin,0); % length of window in the matlab date format

for ii = 1 : numel(windowStart); % number of windows
    
    fprintf('\n Xcor window start time: %s\n',datestr(absStartTimes(ii)));
    
    test = find( abs( absStartTimes(ii) - matlabDate ) <= windowMatlab);    
    
    if isempty(test)
        fprintf('\t No earthquakes near window\n');
    else
        
        % Determine if EQ is before or after window this way we don't throw
        % out windows we do not need to. We will need Piero to fix the
        % CatalogParser2.m so that we know the minutes of the EQs.
        
        %         idx = ( absStartTimes(ii) - matlabDate(test) ) <= 0; % (negative if EQ began after window)
        %         distTime = absStartTimes(ii)-matlabDate(test);

        % For now just plot the EQ information if it occurs within the xcor
        % window
        for jj = 1:numel(test)
                fprintf('\t EQumber %d, EQtime %s, EQsize %2.2f.\n',Catalog.number(test(jj)),char(Catalog.time2(test(jj))),Catalog.magnitude(jj));      
        end
        
        
    end
   
    % Besides throwing out xcor time windows containing EQs, we could also
    % think about setting some kind of magnitude and distance threshold.
    % For example, if a magnitude 4 EQ occurs far from the stations, should
    % we remove those correlations? Or should we only worry about EQs with
    % M > 5? These can be explored more by Celeste.
    
end























