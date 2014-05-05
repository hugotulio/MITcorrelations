clear all
close all
clc

% This script creates day long matrices of all seismic data data in the
% directory USArray

lonDir = '/Volumes/ChainRing/USArrayData/'; % directory containing data subdivided by longitude

lonFolder      = dir(lonDir); % list of data folders
lonFolder(1:2) = []; % get rid of '.' and '..' in dir() output
nFolder        = numel(lonFolder);      % number of folders

% process folders
for ii = 1:nFolder
    fprintf('Checking folder %s\n',lonFolder(ii).name);
    tmp_ii = [lonDir lonFolder(ii).name];
    dateFolder = dir(tmp_ii); % get folder list
    dateFolder(1:2) = []; % get rid of '.' and '..' in dir() output
    
    for  jj = 1:numel(dateFolder);
        fprintf('Checking folder %s\n',dateFolder(jj).name);
        tmp_jj = [tmp_ii '/' dateFolder(jj).name];
        continuousFolder = dir(tmp_jj); % get folder list
        continuousFolder(1:2) = []; % get rid of '.' and '..' in dir() output
        
        for kk = 1:numel(continuousFolder)
            fprintf('Checking folder %s\n',continuousFolder(kk).name);
            tmp_kk = [tmp_jj '/' continuousFolder(kk).name '/BH/'];
            
            sacFile      = dir(tmp_kk); % list of deconvolved SAC files
            sacFile(1:2) = []; % get rid of '.' and '..' in dir() output
            
            for ll = 1:numel(sacFile)
                % load the ObspyDMT deconvolved data
                w = loadsacfile({[tmp_kk sacFile(ll).name]});
                
                plot(w);
                
            end
        end
    end
end


