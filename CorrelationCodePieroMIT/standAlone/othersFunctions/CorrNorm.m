% normalized correlation
% calculation of noise correlation using 4 different normalization
% approaches see Bendat & Piersol 2000
% 
% 1) array normalization (Transfer function) using estimated spectra based on array observation
% 2) station normalization (Transfer function) spectrum must be smoothed
% 3) simple normalization (Coherence)
% 4) Autocorr energy normalization
% 
% INPUT PARAMETERS:
% 
% s1 and s2 : traces to be treated
% Wn: time bandwidth product
% 
% OUTPUT:
% 
% c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
% c2: simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
% c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
% c4: Transfer function array normalization C12(w)/({abs(Sarray(w))^2})
%
% 
% CREATED BY Piero Poli Massachussets Institute of Technology
% V.1.1 - 18 Sept 2013


function [c1,c2,c3]=CorrNorm(s1,s2,Fs,SmoothMethod,Wn,K)


% general formulation

corrLength=length(s1)+length(s2)-1; % get the final correlation length

SpectralDomainCorr=fft(s1,corrLength).*conj(fft(s2,corrLength)); % calculate correlation

% 1) get the normalized correlation as in Bendat Piersol pp. 125 eq. 5.16

sigma1 = fftshift(ifft(fft(s1,corrLength).*conj(fft(s1,corrLength)))); %autocorr s1
sigma2 = fftshift(ifft(fft(s2,corrLength).*conj(fft(s2,corrLength)))); %autocorr s2

norm1=sigma1((length(sigma1)+1)/2).*ones(1,corrLength); % get t0 the enrgy for autocorr s1
norm2=sigma1((length(sigma2)+1)/2).*ones(1,corrLength); % get t0 the enrgy for autocorr s2

c1 = fftshift(ifft(SpectralDomainCorr))./(sqrt(norm1).*sqrt(norm2)); % normalize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      DONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spectral smoothing
if strcmp(SmoothMethod,'taper')==1 %spectral smoothing using taper and sepian

    [e,v] = dpss(length(s1),Wn,K); % compute the slepian sequences 

[ns1,f] = pmtm(s1',e,v,corrLength,Fs,'DropLastTaper',false,'twosided'); % get smoothed spectrum for s1
[ns2,f] = pmtm(s2',e,v,corrLength,Fs,'DropLastTaper',false,'twosided'); % get smoothed spectrum for s2

% get time domain normalized correlation
c2 = fftshift(ifft(SpectralDomainCorr./sqrt(ns1'.*ns2'))).*tukeywin(corrLength,0.1)';

% get time domain normalized correlation by station
c3 = fftshift(ifft(SpectralDomainCorr./(ns1'.^2)));

elseif strcmp(SmoothMethod,'median')==1

n1 = abs(fft(s1,corrLength)); % get abs spectrum for s1   
n2 = abs(fft(s2,corrLength)); % get abs spectrum for s2

ns1 = medfilt1(n1,Wn)';
ns2 = medfilt1(n2,Wn)';

% get time domain normalized correlation
c2 = real(fftshift(ifft(SpectralDomainCorr./(ns1'.*ns2')))).*tukeywin(corrLength,0.1)';

% get time domain normalized correlation by station
c3 = real(fftshift(ifft(SpectralDomainCorr./(ns1'.^2)))).*tukeywin(corrLength,0.1)';

disp('.')

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      DONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%

%c2 = fftshift(ifft(SpectralDomainCorr./(sigma_12.*sigma_22)));
% 1) array normalization using estimated spectra based on array observation


% 2) station normalization (Transfer function) spectrum must be smoothed

% 3) simple autospectra normalization

