clear all
close all
clc

addpath('/Users/dmikessell/GIT/MITcorrelations/matlab');

% This example script shows how to sort the time-window correlation
% matrices from runDayCorrelations.m into matrices for each stations
% pair containing all of the time windows. The newly created matrices
% can then be processed and stacked to estimate the Green's function.

inputDirectory = './4hour_test';
outputDirectory = './4hour_test_sorted_append';

[success, message] = rmdir(outputDirectory,'s');

% sort the correlations by station pairs for all time windows
tic
% sortDayCorrelations( inputDirectory, outputDirectory);
toc

%%


outputDirectory = './4hour_test_sorted_append';

tic
sortDayCorrelations3( inputDirectory, outputDirectory);
toc

%%

matObj = matfile([outputDirectory '/M26A-HWUT.mat'],'Writable',true);

%%  try to write a waveform object by appending to existing file


W = waveform(); % make a single waveform
S = struct(W);
save('wave.mat', 'S', '-v7.3'); % save matfile
matObj = matfile('wave.mat', 'Writable', true); % open matfile and make writable
nWave = numel(matObj.S); % number of waveforms in matfile already
matObj.S(1,nWave + 1) = S; % append another waveform

%%

close all
clear all

load('wave.mat');
