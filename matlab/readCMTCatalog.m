function out = readCMTCatalog(catalog, plot, year)
%
% USAGE: out = readCMTCatalog(catalog, plot, year)
%
% Read seismic catalog fom CMT project 'List of Events Name'
%
% INPUT:
%   catalog = catalog name, including path
%   plot = 'yes' to plot matrix of EQ occurence
%   year = choose a year
% OUTPUT
%   out = a structure with the following fields for each EQ
%       out.year
%       out.month
%       out.day
%       out.hour
%       out.minute
%
% Piero Poli Massachussets Institute of Technology email: ppoli@mit.edu
% 18 Sept 2013
% Last modified by Dylan Mikesell
% 27 June 2014

warning off

[catalog.date, catalog.region] = textread(char(catalog),'%s %s %*[^\n]');

for ii = 1 : size(catalog.date,1)
    out.year(ii,:)   = catalog.date{ii}(1:4);
    out.month(ii,:)  = catalog.date{ii}(5:6);
    out.day(ii,:)    = catalog.date{ii}(7:8);
    out.hour(ii,:)   = catalog.date{ii}(9:10);
    out.minute(ii,:) = catalog.date{ii}(11:12);
end

% get julian date
DATE = [ str2num(out.year) , str2num(out.month) , str2num(out.day) ];

for ii = 1 :size(DATE,1)
    J(ii,:) = julian( DATE(ii,1), DATE(ii,2), DATE(ii,3)) - julian(year,1,1) + 1;
end
out.julianD = J;

%--------------------------------------------------------------------------
% plot matrix data of eq occurrence
if strcmp(plot,'yes') == 1
    
    plot((str2num(out.hour)+1),J,'ro')
    xlabel('hours')
    ylabel('julian day')
end

return