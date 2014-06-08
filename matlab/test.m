clear all
close all
clc

addpath('/Users/dmikessell/GIT/MITcorrelations/matlab');
addpath('/Users/dmikessell/GIT/MITcorrelations/matlab/external');

% load a test matrix
load('/Users/dmikessell/workspace/IRIS/JulianDayData/julDay_1152.mat');
% These data come into MATLAB as waveform objects. This way they
% contain all the META-information we need.

%% setup the correlation output structure

outputDirectory = './output_test';

% make directory for output of the correlations
[success,message,messageID] = checkOutputDir(outputDirectory);
        
%% plot to see raw data

C = correlation(W);
plot(C);

%% remove any blank traces

blankIdx    = ( sum(double(W),1) == 0 ); % find blank traces
W(blankIdx) = []; % remove blank traces

%% frequency filter data

fmin = 0.1; % low-cutt of filter (Hz)
fmax = 1.0; % high-cut of filter (Hz)
btord = 3; % number of poles in Butterworth filter

Filt   = filterobject('B',[fmin,fmax],btord); % create a frequency filter for waveform
W_filt = filtfilt(Filt,W); % apply the filter to the data 

%% remove spikes in data

% This routine divides the time series up in to window and checks for
% spike using an energy test. If the maximum energy in the time
% exceeds a ratio computed by std. deviation of the energy in the
% window, the window, and the precedding window, are set to zero.

WindowMin = 60*4; % window lenght (min)
Nstd = 20; % threshold for determining whether or not a spike.

W_rsp = RemoveSpikesWaveform2(W,WindowMin,Nstd);

%% Preprocess data

% In this section, we need to despike the data and determine whether
% or not we will 'whiten' the data. If we whiten, we need to choose
% the method for whitening.

FB = [0.1, 1];
W_whiten = WhitenWaveform(W,FB);

%% Crosscorrelate all combinations for data.

% We do this by first computing the frequency domain data that we need.

% Next we do a double loop over all possible stations implemented in
% parallel.

C = runDayCorrelations(W,winLength,overlapPercent);


%% Postprocess the correlations

% Here we stack the data using different  methods


