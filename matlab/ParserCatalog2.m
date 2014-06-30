function Catalog = ParserCatalog2(CATALOG)
%
% USAGE: Catalog = ParserCatalog2(CATALOG)
%
%   parser to read the EVENT-INFO CATALOG from obspy-DMT
%   !!!! THIS VERSION LOAD THE CATALOG FROM NEW OBSPY-DMT VERSION !!!!
%
% INPUT:
%   CATALOG = string filename, including path to EVENT-CATALOG
% OUTPUT:
%   Catalog = a structure with the following fields
%       number:    ID of the earthquake
%       time:      date of earthquake
%       depth:     hyponcetral depth
%       lon:       hypocenter lon
%       lat:       hypocenter lat
%       magnitude: magnitude from catalog
%
% CREATED BY Piero Poli Massachussets Institute of Technology
% V.1.2 - 2 April 2014
% Modified by Dylan Mikesell (MIT) 27 June 2014

% Find all attributes
fid1       = fopen(char(CATALOG),'r');
attributes = textscan(fid1,'%*s %s','CommentStyle','%','CollectOutput',1,'Headerlines',17,'delimiter',':','MultipleDelimsAsOne', 1);
attributes = attributes{1};
fclose(fid1);

% Select Attributes
lines = 1 : 11 : size(attributes,1)-9;

for ii = 1 : numel(lines)
    Catalog.number(ii)        = str2double(char(attributes(lines(ii))));
    Catalog.catalogType(ii,:) = (attributes(lines(ii)+1));
    Catalog.magnitude(ii)     = str2double(char(attributes(lines(ii)+5)));
    Catalog.depth(ii)         = str2double(char(attributes(lines(ii)+6)));
    Catalog.lat(ii)           = str2double(char(attributes(lines(ii)+7)));
    Catalog.lon(ii)           = str2double(char(attributes(lines(ii)+8)));
    Catalog.time1(ii)         = (attributes(lines(ii)+2));
    Catalog.time2(ii)         = (attributes(lines(ii)+3));
end

return