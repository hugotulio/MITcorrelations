% normalized correlation with input spectra
% calculation of noise correlation using 4 different normalization
% approaches see Bendat & Piersol 2000
% 
% 1) array normalization (Transfer function) using estimated spectra based on array observation
% 
% INPUT PARAMETERS:
% 
% s1 and s2 : traces to be treated
% Wn: time bandwidth product
% 
% OUTPUT:
% 
% c4: Transfer function array normalization C12(w)/({abs(Sarray(w))^2})
%
% 
% CREATED BY Piero Poli Massachussets Institute of Technology
% V.1.1 - 18 Sept 2013


function [c4,norm]=CorrNormArray(s1,s2,SmoothMethod,spectrum,Wn)

 
% general formulation

corrLength=length(s1)+length(s2)-1; % get the final correlation length

SpectralDomainCorr=fft(s1,corrLength).*conj(fft(s2,corrLength)); % calculate correlation


if strcmp(SmoothMethod,'median')==1


norm = medfilt1(spectrum,Wn)';


% get time domain normalized correlation by station
c4 = real(fftshift(ifft(SpectralDomainCorr./(norm'.^2)))).*tukeywin(corrLength,0.1)';
disp('.')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      DONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
