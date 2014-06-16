clear all
close all
clc

addpath('/Users/dmikessell/GIT/MITcorrelations/matlab');

% This example script shows how to sort the time-window correlation
% matrices from runDayCorrelations.m into matrices for each stations
% pair containing all of the time windows. The newly created matrices
% can then be processed and stacked to estimate the Green's function.

inputDirectory = './4hour_test';
outputDirectory = './4hour_test_sorted';

[success, message] = rmdir(outputDirectory,'s');

% sort the correlations by station pairs for all time windows
sortDayCorrelations( inputDirectory, outputDirectory);


