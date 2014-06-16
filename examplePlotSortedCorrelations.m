clear all
close all
clc

% This script sorts and plots a matrix of station-pair correlations.

inputDirectory = './4hour_test_sorted';

files = dir([inputDirectory '/*.mat']);

idx = 1;
files(idx).name

load([inputDirectory '/' files(idx).name]);

%%

data = cell2mat(get(statC,'c1')); % get c1 data
data = reshape(data,[numel(get(statC(1),'c1')),numel(statC)]); % make c1 vector into matrix

plot(sum(data,2))

% imagesc(data)
% wigb(data)
