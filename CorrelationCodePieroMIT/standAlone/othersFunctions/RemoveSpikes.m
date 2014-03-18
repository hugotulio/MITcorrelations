% function to remove EQs signal from seismic time series.
%
% PROCESSING:
% 1) load single station data
% 2) filtering
% 3) windowing data
% 4) remove eqs using catalog
% 5) calculate energy of seismic noise using no-null remaining window
% 6) evaluate the STD of the no-null daily energy (should be a level of noise energy for that day) 
% 7) chek if window has max energy larger than N times the daily energy STD
% 8) cosine taper the remaining no-null windows
%
% 
% PARAMETERS:
% -Window: size of the window in minute.
%
% -catalog: char value of the EQ catalog used. Function to read catalog is
% at the end of this fucntion. Catlog format is: 200701010031A CARLSBERG
% RIDGE  (eq name GMT output http://www.globalcmt.org/CMTsearch.html)
%
% -Fs: sample frequency
%
% -Nstd: integer for defining std threshold
%
% CREATED BY Piero Poli Massachussets Institute of Technology
% V.1.1 - 19 Sept 2013

function trace = RemoveSpikes(inputTrace,catalog,Fs,Window,Nstd,FB,JulDay)

trace = inputTrace;

% 1) load single station data
trace(end+1)=0; % add a sample to ensure the windowing being right. Sample be removed at the end.

% 2) filtering
[a1,b1]=butter(2,FB*2/Fs);
traceFilt = filtfilt(a1,b1,trace);
% 3) windowing data

    % 3.1 define the windows
    Win = Window * 60 * Fs;
    t0win = 1 : Win : (86400*Fs+1 - Win) ;
 
    % 3.2 Get information from seismic catalog
    out=readCMTCatalog(catalog,'no',2007);
        % 3.2.1 find earthquakes in this day
        A = find(out.julianD==JulDay);
        % 3.2.2 get time of earthquakes
        houreqs=(str2num(out.hour(A,:)))*3600*Fs;
        minuteeqs=(str2num(out.minute(A,:)))*60*Fs;
        EqTime=houreqs+minuteeqs;
        clear A

        % window analysis of eqs presence
        for i2 = 1 : numel(t0win) - 1
            
            indx = find(EqTime>t0win(i2)&EqTime<t0win(i2+1), 1); 
            
            % 4) remove eqs using catalog
            if isempty(indx) == 0
                if i2 > 21
                trace(t0win(i2):end)=0;
                else
                trace(t0win(i2):t0win(i2+2))=0;
                end
            end
        end

% 5) calculate energy of seismic noise using no-null remaining window
  
ee = (traceFilt.^2)>0; % make energy non-null tace
Mu = std(traceFilt(ee).^2);  % get std for non-null trace (daily value of  noise energy...)
       

% 7) chek if window has max energy larger than N times the daily energy STD

 for i3 = 1 : numel(t0win) -1
     
    win=traceFilt(t0win(i3):(t0win(i3)+Win)); 
    
    wabs(i3)=max((win).^2);
    
    if wabs(i3)>Nstd*Mu
    
    if i3 ==1
    trace(t0win(i3):(t0win(i3)+Win))=0;    
    else
    trace(t0win(i3)-Win:(t0win(i3)+Win))=0;
    end
    
    else
        
    
    Tw = tukeywin(Win+1,.1);    
    trace(t0win(i3):(t0win(i3)+Win))=trace(t0win(i3):(t0win(i3)+Win)).* Tw;
    
    % check how many sample in the window are larger than zero
    indx2 = find(abs(trace(t0win(i3):(t0win(i3)+Win)))>1e-40);
    
    % remove windows with zeros
    if numel(indx2)<Win-1000
        trace(t0win(i3):(t0win(i3)+Win))=0;   
    end

    end
    clear win 
 end
        
trace(end)=[]; % remove the sample added at the beginning...

%%%%%% OTHER USEFUL FUNCTIONS....

% read seismic catalog fom CMT project 'List of Events Name'
% Piero Poli Massachussets Institute of Technology email: ppoli@mit.edu
% 18 Sept 2013

function [out]=readCMTCatalog(catalog,plot,year)
warning off
[catalog.date, catalog.region]=textread(char(catalog),'%s %s %*[^\n]');

for i1 = 1 :  size(catalog.date,1)

out.year(i1,:) = catalog.date{i1}(1:4);
out.month(i1,:) = catalog.date{i1}(5:6);
out.day(i1,:) = catalog.date{i1}(7:8);
out.hour(i1,:) = catalog.date{i1}(9:10);
out.minute(i1,:) = catalog.date{i1}(11:12);
end


% get julian date
DATE = [str2num(out.year) , str2num(out.month) , str2num(out.day)];


for i2 = 1 :size(DATE,1)
J(i2,:) = julian(DATE(i2,1),DATE(i2,2),DATE(i2,3))-julian(year,1,1)+1;
end
out.julianD = J;

% plot matrix data of eq occurrence

if strcmp(plot,'yes')==1

plot((str2num(out.hour)+1),J,'ro')
xlabel('hours')
ylabel('julian day')
end

function c = julian(yy,mm,dd,hh,m,s)

% JULIAN  Converts Gregorian calendar dates to Julian day number
%
% Although the formal definition holds that Julian days start 
%  and end at noon, here Julian days start and end at midnight.
%
% In this convention, Julian day 2440000 began at 0000 hours, May 23, 1968.
%
%  Converts Date to decimal Julian Day:
%
%    3 .. 6 Inputs:
%
%      JulianDay = JULIAN( YY , MM , DD , [hh] , [mm] , [ss] )
%
%    single Input, 3 .. 6 Columns:
%
%      JulianDay = JULIAN([ YY MM DD [hh] [mm] [ss] ]);
%
%  Converts Julian Day to Date, single Input with 1 Column:
%
%    [YY MM DD hh mm ss] = JULIAN( JulianDay )
%
%     ************************************************************
%
%        YY.... year (e.g., 1979) component
%        MM.... month (1-12) component
%        DD.... day (1-31) component of Gregorian date
%        hh.... decimal hours (assumed 0 if absent)
%        mm.... minutes
%        ss.... seconds
%
%        JulianDay  decimal Julian Day
%
%     ************************************************************
%
% see also: DATENUM, DATEVEC
%
%     recoded for MATLAB  by Rich Signell, 5-15-91
%     improved by CBegler 02/2005
%

Nin = nargin;

siz = [ 1 3 4 5 6 ];

if ~any( Nin == siz )
    error('Invalid Number of Inputs.');
end

if Nin == 1

   if isempty(yy)
      c = [];
      return
   end

   s2 = size(yy,2);

   if ~any( s2 == siz )
       error('Invalid Size of Input.');
   end

   %------------------------------------------------------
   % Day --> Date
   %------------------------------------------------------
   if s2 == 1
      c = datetime( yy + datenum(1968,05,23) - 2440000 );
      return
   end
   %------------------------------------------------------

   yy = cat( 2 , y , zeros(size(yy,1),6-size(yy,2)) );

   hh = yy(:,4) + yy(:,5)/60 + yy(:,6)/3600;
   dd = yy(:,3);
   mm = yy(:,2);
   yy = yy(:,1);

else

   if Nin < 4, hh = 0; end 
   if Nin < 5, m  = 0; end 
   if Nin < 6, s  = 0; end 

   hh = hh + m/60 + s/3600;

   sy = size(yy); py = prod(sy);
   sm = size(mm); pm = prod(sm);
   sd = size(dd); pd = prod(sd);
   sh = size(hh); ph = prod(sh);

   if ~( ( isequal(sy,sm) | ( py == 1 ) | ( pm == 1 ) ) & ...
         ( isequal(sy,sd) | ( py == 1 ) | ( pd == 1 ) ) & ...
         ( isequal(sy,sh) | ( py == 1 ) | ( ph == 1 ) ) & ...
         ( isequal(sm,sd) | ( pm == 1 ) | ( pd == 1 ) ) & ...
         ( isequal(sm,sh) | ( pm == 1 ) | ( ph == 1 ) ) & ...
         ( isequal(sd,sh) | ( pd == 1 ) | ( ph == 1 ) )       )
       error('Matrix Dimensions must be agree.');
   end

end

mo=mm+9;
yr=yy-1;

ii = ( mm > 2 );
if any(ii)
       ii  = find(ii);
    mo(ii) = mm(ii) - 3;
    yr(ii) = yy(ii);
end
 
c = floor(yr/100);

yr = yr - c*100;
      
c = floor((146097*c)/4) + floor((1461*yr)/4) + ...
    floor( ( 153*mo + 2 ) / 5 ) + dd + 1721119;

%     If you want julian days to start and end at noon, 
%     replace the following line with:
%     j=j+(h-12)/24;
 
c = c + hh/24;

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function d = datetime(t);

% DATETIME  returns correct DateTime
%
% Takes care on accuraccy-problems of DATEVEC, 
%    which returns seconds == 59.99999999999272
%

d      = datevec(t);

d(:,6) = round(d(:,6));  % Round seconds

dd = d(:,3);  % Original DayNumber

quot = [ 60 60 24 ]; %  [ ss-->mm  mm-->hh  hh-->dd ]
ind  = [ 6  5  4  ];

for ii = 1 : 3

    p = fix( d(:,ind(ii)) / quot(ii) );

 d(:,ind(ii)-0) = d(:,ind(ii)-0) - p * quot(ii);
 d(:,ind(ii)-1) = d(:,ind(ii)-1) + p;
  
end

% Check if DayNumber has changed

ii = find( d(:,3) > dd );

if isempty(ii)
   d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];
   return
end

% New Date

[d(ii,1),d(ii,2),d(ii,3)] = datevec( datenum(d(ii,1),d(ii,2),d(ii,3)) );


d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];

