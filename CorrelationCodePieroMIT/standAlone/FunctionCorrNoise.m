% noise cross-correlation code

% CREATED BY Piero Poli Massachussets Institute of Technology
% V.1.1 - 19 Sept 2013

addpath /Users/pieropoli/NoiseCorrelation/CorrelationCode/standAlone/othersFunctions

%function FunctionCorrNoise()

disp('SCC (Seismic Correlation Code) - CREATED BY Piero Poli Massachussets Institute of Technology')

%%%%%%%%%%%%%%  PARAMETERS

tic
InpDir= input('Input directoy of data: ','s');

tic
OutDir= input('Output directoy of data: ','s');

tic
stationInfolist = input('Reference Stations List : ','s');

tic
arrayInfolist = input('Array Stations List : ','s');

tic
SpectralSmoothing = input('Spectral Smoothing Type : ','s');

tic
SpectralSmoothingValue = input('Spectral Smoothing Value : ');

tic
SpectralSmoothingValueArray = input('Spectral Smoothing Value Array : ');

tic
minFB_EQ = input('minimum freq. Filtering FB for remove earthquakes : ');

tic
maxFB_EQ = input('maximum freq. Filtering FB for remove earthquakes : ');

tic
minFB = input('minimum freq. Filtering FB : ');

tic
maxFB = input('maximum freq. Filtering FB : ');

tic
Nstd = input('std limit to remove EQS : ');

tic
Window = input('Size of windows to calculate noise correlations (in minutes): ');

tic
Fs = input('Sample frequency (hz) : ');

tic
year = input('Enter the year : ');

tic
catalog = input('earthquake catalog (GMT) : ','s');

% OK SETTING PARAMETERS

FB_EQ = [minFB_EQ maxFB_EQ];

FB = [minFB maxFB ];

mkdir(OutDir)

% 1) read the goegraphical info of the stations

[stationInfo.net,stationInfo.sta,stationInfo.lat,stationInfo.lon]=textread(char(stationInfolist),'%s %s %f %f %*[^\n]');

[arrayInfo.net,arrayInfo.sta,arrayInfo.lat,arrayInfo.lon]=textread(char(arrayInfolist),'%s %s %f %f %*[^\n]');


% 3) Load and process data to remove spikes


% may preallocate a matrix with dayily race to be saved ....

day = 1 : 365;

    for i2 = day % loop over day
    
        
        
        
        if size(int2str(i2),2)==1
        DAY=char(['00' int2str(i2)]);    
        elseif size(int2str(i2),2)==2
        DAY=char(['0' int2str(i2)]);    
        else
        DAY=int2str(i2);
        end
        
        % make directory for output of the correlations
        eval(['mkdir ' char(OutDir) '/' char(DAY)])
        
        %%%% now load station data...
        % preallocate s1
        
        time1 = cputime;
        for i1 = 1 : size(stationInfo.net,1) % load the station file
        
            if exist([char(InpDir) '/' DAY '/' char(stationInfo.sta(i1,:)) '.mat'])==2
            disp([char(InpDir) '/' DAY '/' char(stationInfo.sta(i1,:)) '.mat'])
            load([char(InpDir) '/' DAY '/' char(stationInfo.sta(i1,:)) '.mat'])
            % treat data for the station
            
            trace1 = trace;clear trace
            char(catalog)
            Fs
            Window
            Nstd
            FB_EQ
            i2
            trace1
            s1(i1,:) = RemoveSpikes(trace1,char(catalog),Fs,Window,Nstd,FB_EQ,i2);
            
            % make the header for stations
            headerStation.net(i1,:) = stationInfo.net(i1,:);
            headerStation.lat(i1,:) = stationInfo.lat(i1,:);
            headerStation.lon(i1,:) = stationInfo.lon(i1,:);
            headerStation.sta(i1,:) = stationInfo.sta(i1,:);
            end
        end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%% now load array data...
        %preacllocate s2
        
        
        for i3 = 1 : size(arrayInfo.net,1) % load the station file
        
            if exist([InpDir '/' DAY '/' char(arrayInfo.sta(i3,:)) '.mat'])==2

            load([InpDir '/' DAY '/' char(arrayInfo.sta(i3,:)) '.mat']);
            % treat data for the array stations
            s2(i3,:) = RemoveSpikes(trace,char(catalog),Fs,Window,Nstd,FB_EQ,i2);
            % make the header for stations
            headerArray.net(i3,:) = arrayInfo.net(i3,:);
            headerArray.lat(i3,:) = arrayInfo.lat(i3,:);
            headerArray.lon(i3,:) = arrayInfo.lon(i3,:);
            headerArray.sta(i3,:) = arrayInfo.sta(i3,:);
            end
        
        end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

% check that both s1 and s2 contain signals
       
       if exist('s1')==0||exist('s2')==0
          
          disp('no data') 
          
       else  
            
        % remove null data from s1
        
        check = mean(s1,2)==0;
        s1(check,:)=[];
        
        % remove header null-trace items
        ns1 = fieldnames(headerStation);
        
        for hdr1 = 1 : size(ns1,1)
            
           eval(['headerStation.' char(ns1(hdr1,:)) '(check,:)=[];']) 
            
        end
        
        clear check ns1 hdr1
        
        % remove null data from s2
        
        check = mean(s2,2)==0;
        s2(check,:)=[];
        
        
        % remove header null-trace items
        ns1 = fieldnames(headerArray);
        
        for hdr2 = 1 : size(ns1,1)
            
           eval(['headerArray.' char(ns1(hdr2,:)) '(check)=[];']) 
            
        end
        
        clear check ns1 hdr2
 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
           %%% windowing correlation %%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

        %%% making window estimate spectra for array
        
        %%% making window estimate spectra for array
        Win = Window * 60 * Fs;
        t0win = 1 : Win : (86400*Fs+1 - Win) ;
        corrLength=Win+Win+1;
        AvSpec = zeros(size(t0win,2),corrLength);
        
                for ws1 = 1 : numel(t0win) -1
                    sp=zeros(size(s2,1),corrLength);
                    for cc12 = 1 : size(s2,1)
                    s2win=s2(cc12,t0win(ws1):(t0win(ws1)+Win));
                    sp(cc12,:)=abs(fft(s2win,corrLength)); 
                    idxchk(cc12) = sum(sp(cc12,2));
                    end
                if size(sp,1)==1    
                AvSpec(ws1,:)=sp;clear sp
                else
                ff =find(idxchk==0);clear idxchk;
                sp(ff,:)=[];clear ff
                AvSpec(ws1,:)=mean(sp,1);clear sp ff
                end
                end

        %%%% now make the windowing and correlation
               
        for cc1 = 1 : size(s1,1)
          
           for cc2 = 1 : size(s2,1) 


           Win = Window * 60 * Fs;
           t0win = 1 : Win : (86400*Fs+1 - Win) ;
           
           % preallocate matrixes
           
           c1 = zeros(length(t0win),Win*2+1);c2=c1;c3=c1;c4=c1;
           
           
           for w1 = 1 : numel(t0win) -1
     
             s1win=s1(cc1,t0win(w1):(t0win(w1)+Win)); 
             s2win=s2(cc2,t0win(w1):(t0win(w1)+Win));
             
             
             
             if sum(s1win)~=0&&sum(s2win)~=0
             % correlation
             [c1(w1,:),c2(w1,:),c3(w1,:)]=CorrNorm(s1win,s2win,Fs,char(SpectralSmoothing),SpectralSmoothingValue,0);  
             [c4(w1,:),spec]=CorrNormArray(s1win,s2win,char(SpectralSmoothing),AvSpec(w1,:),SpectralSmoothingValueArray); 
             else 
                disp('no data for that window') 
             end
             clear s1win s2win   
           end
           
           % check the correlations
           check = find(mean(c1,2)==0);
           c1(check,:)=[];
           c2(check,:)=[];
           c3(check,:)=[];
           c4(check,:)=[];
           % average daily data
           C1 = mean(c1,1);
           C2 = mean(c2,1);
           C3 = mean(c3,1);
           C4 = mean(c4,1);
           clear c1 c2 c3 c4
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           
           % make correlation header
           corr.header.st1=headerStation.sta(cc1,:);
           corr.header.st2=headerArray.sta(cc2,:);
           corr.header.lat_st1=headerStation.lat(cc1,:);
           corr.header.lon_st1=headerStation.lon(cc1,:);
           corr.header.lat_st2=headerArray.lat(cc2,:);
           corr.header.lon_st2=headerArray.lon(cc2,:);
           corr.C1=C1;
           corr.C2=C2;
           corr.C3=C3;
           corr.C4=C4;
           corr.header.date=date;
           corr.header.Jday=DAY;
           corr.header.year=year;
           corr.header.InpDir=char(InpDir);
           corr.header.OutDir=char(OutDir);
           corr.header.year=year;
           corr.header.FB=FB;
           corr.header.Window=Window;
           corr.header.SpectralSmoothing=char(SpectralSmoothing);
           corr.header.SpectralSmoothingValue=SpectralSmoothingValue;
           corr.header.SpectralSmoothingValueArray=SpectralSmoothingValueArray;
           corr.header.FB_EQ=FB_EQ;
           corr.header.Nstd=Nstd;
                      
           clear C1 C2 C3 C4
           % here save the correlation for the day and for the array-station
           % combination
           disp(['corr ' char(headerStation.sta(cc1)) ' ' char(headerArray.sta(cc2,:)) ])
           eval(['save ' char(OutDir) '/' char(DAY) '/corr_' char(headerStation.sta(cc1,:)) '_' char(headerArray.sta(cc2,:)) ' corr'])
           clear corr
           
           end

 
        end


       end
        
        %%%%%%%%%%%%%% DONE %%%%%%%%%%%%%%%%%%%
        
        disp(['end of day ' int2str(i2)])
        clear s1 s2
    end
%     
%     
%     
%     
%     %%%% ADDITIONAL FUNCTIONS
%     
% %     function to remove EQs signal from seismic time series.
% % 
% % PROCESSING:
% % 1) load single station data
% % 2) filtering
% % 3) windowing data
% % 4) remove eqs using catalog
% % 5) calculate energy of seismic noise using no-null remaining window
% % 6) evaluate the STD of the no-null daily energy (should be a level of noise energy for that day) 
% % 7) chek if window has max energy larger than N times the daily energy STD
% % 8) cosine taper the remaining no-null windows
% % 
% % 
% % PARAMETERS:
% % -Window: size of the window in minute.
% % 
% % -catalog: char value of the EQ catalog used. Function to read catalog is
% % at the end of this fucntion. Catlog format is: 200701010031A CARLSBERG
% % RIDGE  (eq name GMT output http://www.globalcmt.org/CMTsearch.html)
% % 
% % -Fs: sample frequency
% % 
% % -Nstd: integer for defining std threshold
% % 
% % CREATED BY Piero Poli Massachussets Institute of Technology
% % V.1.1 - 19 Sept 2013
% 
% function trace = RemoveSpikes(inputTrace,catalog,Fs,Window,Nstd,FB,JulDay)
% 
% trace = inputTrace;
% 
% % 1) load single station data
% trace(end+1)=0; % add a sample to ensure the windowing being right. Sample be removed at the end.
% 
% % 2) filtering
% [a1,b1]=butter(2,FB*2/Fs);
% traceFilt = filtfilt(a1,b1,trace);
% % 3) windowing data
% 
%     % 3.1 define the windows
%     Win = Window * 60 * Fs;
%     t0win = 1 : Win : (86400*Fs+1 - Win) ;
%  
%     % 3.2 Get information from seismic catalog
%     out=readCMTCatalog(catalog,'no',2007);
%         % 3.2.1 find earthquakes in this day
%         A = find(out.julianD==JulDay);
%         % 3.2.2 get time of earthquakes
%         houreqs=(str2num(out.hour(A,:)))*3600*Fs;
%         minuteeqs=(str2num(out.minute(A,:)))*60*Fs;
%         EqTime=houreqs+minuteeqs;
%         clear A
% 
%         % window analysis of eqs presence
%         for i2 = 1 : numel(t0win) - 1
%             
%             indx = find(EqTime>t0win(i2)&EqTime<t0win(i2+1), 1); 
%             
%             % 4) remove eqs using catalog
%             if isempty(indx) == 0
%                 if i2 > 21
%                 trace(t0win(i2):end)=0;
%                 else
%                 trace(t0win(i2):t0win(i2+2))=0;
%                 end
%             end
%         end
% 
% % 5) calculate energy of seismic noise using no-null remaining window
%   
% ee = (traceFilt.^2)>0; % make energy non-null tace
% Mu = std(traceFilt(ee).^2);  % get std for non-null trace (daily value of  noise energy...)
%        
% 
% % 7) chek if window has max energy larger than N times the daily energy STD
% 
%  for i3 = 1 : numel(t0win) -1
%      
%     win=traceFilt(t0win(i3):(t0win(i3)+Win)); 
%     
%     wabs(i3)=max((win).^2);
%     
%     if wabs(i3)>Nstd*Mu
%     
%     if i3 ==1
%     trace(t0win(i3):(t0win(i3)+Win))=0;    
%     else
%     trace(t0win(i3)-Win:(t0win(i3)+Win))=0;
%     end
%     
%     else
%         
%     
%     Tw = tukeywin(Win+1,.1);    
%     trace(t0win(i3):(t0win(i3)+Win))=trace(t0win(i3):(t0win(i3)+Win)).* Tw;
%     
%     % check how many sample in the window are larger than zero
%     indx2 = find(abs(trace(t0win(i3):(t0win(i3)+Win)))>1e-40);
%     
%     % remove windows with zeros
%     if numel(indx2)<Win-1000
%         trace(t0win(i3):(t0win(i3)+Win))=0;   
%     end
% 
%     end
%     clear win 
%  end
%         
% trace(end)=[]; % remove the sample added at the beginning...
% 
% %%%%%% OTHER USEFUL FUNCTIONS....
% 
% % read seismic catalog fom CMT project 'List of Events Name'
% % Piero Poli Massachussets Institute of Technology email: ppoli@mit.edu
% % 18 Sept 2013
% 
% function [out]=readCMTCatalog(catalog,plot,year)
% warning off
% [catalog.date, catalog.region]=textread(char(catalog),'%s %s %*[^\n]');
% 
% for i1 = 1 :  size(catalog.date,1)
% 
% out.year(i1,:) = catalog.date{i1}(1:4);
% out.month(i1,:) = catalog.date{i1}(5:6);
% out.day(i1,:) = catalog.date{i1}(7:8);
% out.hour(i1,:) = catalog.date{i1}(9:10);
% out.minute(i1,:) = catalog.date{i1}(11:12);
% end
% 
% 
% % get julian date
% DATE = [str2num(out.year) , str2num(out.month) , str2num(out.day)];
% 
% 
% for i2 = 1 :size(DATE,1)
% J(i2,:) = julian(DATE(i2,1),DATE(i2,2),DATE(i2,3))-julian(year,1,1)+1;
% end
% out.julianD = J;
% 
% % plot matrix data of eq occurrence
% 
% if strcmp(plot,'yes')==1
% 
% plot((str2num(out.hour)+1),J,'ro')
% xlabel('hours')
% ylabel('julian day')
% end
% 
% function c = julian(yy,mm,dd,hh,m,s)
% 
% % JULIAN  Converts Gregorian calendar dates to Julian day number
% %
% % Although the formal definition holds that Julian days start 
% %  and end at noon, here Julian days start and end at midnight.
% %
% % In this convention, Julian day 2440000 began at 0000 hours, May 23, 1968.
% %
% %  Converts Date to decimal Julian Day:
% %
% %    3 .. 6 Inputs:
% %
% %      JulianDay = JULIAN( YY , MM , DD , [hh] , [mm] , [ss] )
% %
% %    single Input, 3 .. 6 Columns:
% %
% %      JulianDay = JULIAN([ YY MM DD [hh] [mm] [ss] ]);
% %
% %  Converts Julian Day to Date, single Input with 1 Column:
% %
% %    [YY MM DD hh mm ss] = JULIAN( JulianDay )
% %
% %     ************************************************************
% %
% %        YY.... year (e.g., 1979) component
% %        MM.... month (1-12) component
% %        DD.... day (1-31) component of Gregorian date
% %        hh.... decimal hours (assumed 0 if absent)
% %        mm.... minutes
% %        ss.... seconds
% %
% %        JulianDay  decimal Julian Day
% %
% %     ************************************************************
% %
% % see also: DATENUM, DATEVEC
% %
% %     recoded for MATLAB  by Rich Signell, 5-15-91
% %     improved by CBegler 02/2005
% %
% 
% Nin = nargin;
% 
% siz = [ 1 3 4 5 6 ];
% 
% if ~any( Nin == siz )
%     error('Invalid Number of Inputs.');
% end
% 
% if Nin == 1
% 
%    if isempty(yy)
%       c = [];
%       return
%    end
% 
%    s2 = size(yy,2);
% 
%    if ~any( s2 == siz )
%        error('Invalid Size of Input.');
%    end
% 
%    %------------------------------------------------------
%    % Day --> Date
%    %------------------------------------------------------
%    if s2 == 1
%       c = datetime( yy + datenum(1968,05,23) - 2440000 );
%       return
%    end
%    %------------------------------------------------------
% 
%    yy = cat( 2 , y , zeros(size(yy,1),6-size(yy,2)) );
% 
%    hh = yy(:,4) + yy(:,5)/60 + yy(:,6)/3600;
%    dd = yy(:,3);
%    mm = yy(:,2);
%    yy = yy(:,1);
% 
% else
% 
%    if Nin < 4, hh = 0; end 
%    if Nin < 5, m  = 0; end 
%    if Nin < 6, s  = 0; end 
% 
%    hh = hh + m/60 + s/3600;
% 
%    sy = size(yy); py = prod(sy);
%    sm = size(mm); pm = prod(sm);
%    sd = size(dd); pd = prod(sd);
%    sh = size(hh); ph = prod(sh);
% 
%    if ~( ( isequal(sy,sm) | ( py == 1 ) | ( pm == 1 ) ) & ...
%          ( isequal(sy,sd) | ( py == 1 ) | ( pd == 1 ) ) & ...
%          ( isequal(sy,sh) | ( py == 1 ) | ( ph == 1 ) ) & ...
%          ( isequal(sm,sd) | ( pm == 1 ) | ( pd == 1 ) ) & ...
%          ( isequal(sm,sh) | ( pm == 1 ) | ( ph == 1 ) ) & ...
%          ( isequal(sd,sh) | ( pd == 1 ) | ( ph == 1 ) )       )
%        error('Matrix Dimensions must be agree.');
%    end
% 
% end
% 
% mo=mm+9;
% yr=yy-1;
% 
% ii = ( mm > 2 );
% if any(ii)
%        ii  = find(ii);
%     mo(ii) = mm(ii) - 3;
%     yr(ii) = yy(ii);
% end
%  
% c = floor(yr/100);
% 
% yr = yr - c*100;
%       
% c = floor((146097*c)/4) + floor((1461*yr)/4) + ...
%     floor( ( 153*mo + 2 ) / 5 ) + dd + 1721119;
% 
% %     If you want julian days to start and end at noon, 
% %     replace the following line with:
% %     j=j+(h-12)/24;
%  
% c = c + hh/24;
% 
% %********************************************************************
% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% 
% function d = datetime(t)
% 
% % DATETIME  returns correct DateTime
% %
% % Takes care on accuraccy-problems of DATEVEC, 
% %    which returns seconds == 59.99999999999272
% %
% 
% d      = datevec(t);
% 
% d(:,6) = round(d(:,6));  % Round seconds
% 
% dd = d(:,3);  % Original DayNumber
% 
% quot = [ 60 60 24 ]; %  [ ss-->mm  mm-->hh  hh-->dd ]
% ind  = [ 6  5  4  ];
% 
% for ii = 1 : 3
% 
%     p = fix( d(:,ind(ii)) / quot(ii) );
% 
%  d(:,ind(ii)-0) = d(:,ind(ii)-0) - p * quot(ii);
%  d(:,ind(ii)-1) = d(:,ind(ii)-1) + p;
%   
% end
% 
% % Check if DayNumber has changed
% 
% ii = find( d(:,3) > dd );
% 
% if isempty(ii)
%    d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];
%    return
% end
% 
% % New Date
% 
% [d(ii,1),d(ii,2),d(ii,3)] = datevec( datenum(d(ii,1),d(ii,2),d(ii,3)) );
% 
% 
% d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];
% 
% 
% 
% % normalized correlation
% % calculation of noise correlation using 4 different normalization
% % approaches see Bendat & Piersol 2000
% % 
% % 1) array normalization (Transfer function) using estimated spectra based on array observation
% % 2) station normalization (Transfer function) spectrum must be smoothed
% % 3) simple normalization (Coherence)
% % 4) Autocorr energy normalization
% % 
% % INPUT PARAMETERS:
% % 
% % s1 and s2 : traces to be treated
% % Wn: time bandwidth product
% % 
% % OUTPUT:
% % 
% % c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
% % c2: simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
% % c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
% % c4: Transfer function array normalization C12(w)/({abs(Sarray(w))^2})
% %
% % 
% % CREATED BY Piero Poli Massachussets Institute of Technology
% % V.1.1 - 18 Sept 2013
% 
% 
% function [c1,c2,c3]=CorrNorm(s1,s2,Fs,SmoothMethod,Wn,K)
% 
% 
% % general formulation
% 
% corrLength=length(s1)+length(s2)-1; % get the final correlation length
% 
% SpectralDomainCorr=fft(s1,corrLength).*conj(fft(s2,corrLength)); % calculate correlation
% 
% % 1) get the normalized correlation as in Bendat Piersol pp. 125 eq. 5.16
% 
% sigma1 = fftshift(ifft(fft(s1,corrLength).*conj(fft(s1,corrLength)))); %autocorr s1
% sigma2 = fftshift(ifft(fft(s2,corrLength).*conj(fft(s2,corrLength)))); %autocorr s2
% 
% norm1=sigma1((length(sigma1)+1)/2).*ones(1,corrLength); % get t0 the enrgy for autocorr s1
% norm2=sigma1((length(sigma2)+1)/2).*ones(1,corrLength); % get t0 the enrgy for autocorr s2
% 
% c1 = fftshift(ifft(SpectralDomainCorr))./(sqrt(norm1).*sqrt(norm2)); % normalize
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%      DONE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % spectral smoothing
% if strcmp(SmoothMethod,'taper')==1 %spectral smoothing using taper and sepian
% 
%     [e,v] = dpss(length(s1),Wn,K); % compute the slepian sequences 
% 
% [ns1,f] = pmtm(s1',e,v,corrLength,Fs,'DropLastTaper',false,'twosided'); % get smoothed spectrum for s1
% [ns2,f] = pmtm(s2',e,v,corrLength,Fs,'DropLastTaper',false,'twosided'); % get smoothed spectrum for s2
% 
% % get time domain normalized correlation
% c2 = fftshift(ifft(SpectralDomainCorr./sqrt(ns1'.*ns2'))).*tukeywin(corrLength,0.1)';
% 
% % get time domain normalized correlation by station
% c3 = fftshift(ifft(SpectralDomainCorr./(ns1'.^2)));
% 
% elseif strcmp(SmoothMethod,'median')==1
% 
% n1 = abs(fft(s1,corrLength)); % get abs spectrum for s1   
% n2 = abs(fft(s2,corrLength)); % get abs spectrum for s2
% 
% ns1 = medfilt1(n1,Wn)';
% ns2 = medfilt1(n2,Wn)';
% 
% % get time domain normalized correlation
% c2 = real(fftshift(ifft(SpectralDomainCorr./(ns1'.*ns2')))).*tukeywin(corrLength,0.1)';
% 
% % get time domain normalized correlation by station
% c3 = real(fftshift(ifft(SpectralDomainCorr./(ns1'.^2)))).*tukeywin(corrLength,0.1)';
% 
% disp('.')
% 
% end
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%      DONE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%
% 
% %c2 = fftshift(ifft(SpectralDomainCorr./(sigma_12.*sigma_22)));
% % 1) array normalization using estimated spectra based on array observation
% 
% 
% % 2) station normalization (Transfer function) spectrum must be smoothed
% 
% % 3) simple autospectra normalization
% 
% % normalized correlation with input spectra
% % calculation of noise correlation using 4 different normalization
% % approaches see Bendat & Piersol 2000
% % 
% % 1) array normalization (Transfer function) using estimated spectra based on array observation
% % 
% % INPUT PARAMETERS:
% % 
% % s1 and s2 : traces to be treated
% % Wn: time bandwidth product
% % 
% % OUTPUT:
% % 
% % c4: Transfer function array normalization C12(w)/({abs(Sarray(w))^2})
% %
% % 
% % CREATED BY Piero Poli Massachussets Institute of Technology
% % V.1.1 - 18 Sept 2013
% 
% 
% function [c4,norm]=CorrNormArray(s1,s2,SmoothMethod,spectrum,Wn)
% 
%  
% % general formulation
% 
% corrLength=length(s1)+length(s2)-1; % get the final correlation length
% 
% SpectralDomainCorr=fft(s1,corrLength).*conj(fft(s2,corrLength)); % calculate correlation
% 
% 
% if strcmp(SmoothMethod,'median')==1
% 
% 
% norm = medfilt1(spectrum,Wn)';
% 
% 
% % get time domain normalized correlation by station
% c4 = real(fftshift(ifft(SpectralDomainCorr./(norm'.^2)))).*tukeywin(corrLength,0.1)';
% disp('.')
% 
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%      DONE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%
% 
% 
% 
% 
% 
