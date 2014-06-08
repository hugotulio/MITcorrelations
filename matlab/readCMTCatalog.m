function [out]=readCMTCatalog(catalog,plot,year)
%
% read seismic catalog fom CMT project 'List of Events Name'
% Piero Poli Massachussets Institute of Technology email: ppoli@mit.edu
% 18 Sept 2013
%

warning off

[catalog.date, catalog.region] = textread(char(catalog),'%s %s %*[^\n]');

for i1 = 1 : size(catalog.date,1)
    
    out.year(i1,:)   = catalog.date{i1}(1:4);
    out.month(i1,:)  = catalog.date{i1}(5:6);
    out.day(i1,:)    = catalog.date{i1}(7:8);
    out.hour(i1,:)   = catalog.date{i1}(9:10);
    out.minute(i1,:) = catalog.date{i1}(11:12);
end

% get julian date
DATE = [ str2num(out.year) , str2num(out.month) , str2num(out.day) ];

for i2 = 1 :size(DATE,1)
    J(i2,:) = julian( DATE(i2,1), DATE(i2,2), DATE(i2,3)) - julian(year,1,1) + 1;
end
out.julianD = J;

% plot matrix data of eq occurrence

if strcmp(plot,'yes') == 1
    
    plot((str2num(out.hour)+1),J,'ro')
    xlabel('hours')
    ylabel('julian day')
end

return