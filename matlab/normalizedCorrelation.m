function [c1,c2,c3] = normalizedCorrelation(s1, s2, Fs, SmoothMethod, Wn, K)
%
% This function computes crosscorrelation of s1 and s2. The returned
% data each have a different normalization applied. The computations
% use 3 different approaches. See Bendat & Piersol 2000
%
% USAGE: [c1,c2,c3] = normalizedCorrelation(s1, s2, Fs, SmoothMethod, Wn, K)
%
% INPUT:
%   s1 = first trace (i.e., the virtual source)
%   s2 = second trace
%   Fs = sample frequency (Hz)
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
% Modified by Dylan Mikesell (mikesell@mit.edu)
% Last modified 2 June 2014

% general formulation

npts = numel(s1);

nfft = 2*npts - 1; % get the final correlation length

% Do the FFT once up front for each trace
s1FFT = fft(s1, nfft);
s2FFT = fft(s2, nfft);

%---------------------------------------------------------------------
% c1: Autocorr energy normalized correlation: C12(t)/(C11(0)C22(0))
%---------------------------------------------------------------------

SpectralDomainCorr = s1FFT.*conj(s2FFT); % correlation in frequency domain

% 1) get the normalized correlation as in Bendat Piersol pp. 125 eq. 5.16
sigma1 = fftshift( ifft( s1FFT.*conj(s1FFT) ) ); % autocorr s1
sigma2 = fftshift( ifft( s2FFT.*conj(s2FFT) ) ); % autocorr s2

norm1=sigma1( (length(sigma1) + 1) / 2 ).*ones(1, nfft); % get t=0 enrgy for autocorr s1
norm2=sigma1( (length(sigma2) + 1) / 2 ).*ones(1, nfft); % get t=0 enrgy for autocorr s2

c1 = fftshift( ifft(SpectralDomainCorr) )./( sqrt(norm1).*sqrt(norm2) ); % normalize

%---------------------------------------------------------------------
% c2: Simple normalization (Coherence) C12(w)/({abs(S1(w))}{abs(S2(w))})
% and
% c3: Transfer function station normalization C12(w)/({abs(S1(w))^2})
%---------------------------------------------------------------------

% spectral smoothing
switch SmoothMethod

    case 'taper'
        
        [e,v] = dpss(length(s1), Wn, K); % compute the slepian sequences
        
        ns1 = pmtm(s1', e, v, nfft, Fs, 'DropLastTaper', false, 'twosided'); % get smoothed spectrum for s1
        ns2 = pmtm(s2', e, v, nfft, Fs, 'DropLastTaper', false, 'twosided'); % get smoothed spectrum for s2
        
        % get time domain normalized correlation
        c2 = fftshift( ifft( SpectralDomainCorr./sqrt(ns1'.*ns2') ) ).*tukeywin(nfft,0.1)';
        
        % get time domain normalized correlation by virtual source station
        c3 = fftshift( ifft( SpectralDomainCorr./(ns1'.^2) ) );
        
    case 'median' % apply median smoothing to the freq domain amplitudes before normalization
    
        n1 = abs(s1FFT); % get abs spectrum for s1
        n2 = abs(s2FFT); % get abs spectrum for s2
        
        ns1 = medfilt1(n1, Wn)';
        ns2 = medfilt1(n2, Wn)';
        
        % get time domain normalized correlation
        c2 = real( fftshift( ifft( SpectralDomainCorr./(ns1'.*ns2') ) ) ).*tukeywin(nfft,0.1)';
        
        % get time domain normalized correlation by virtual source station
        c3 = real( fftshift( ifft( SpectralDomainCorr./(ns1'.^2) ) ) ).*tukeywin(nfft,0.1)';
  
end

return