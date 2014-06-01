function [aMatF,fArray] = fft_traces(aMat,dt,fmin,fmax,df)
%
% USAGE: [aMatF,fArray] = fft_traces(aMat,dt,fmin,fmax,df)
%
% INPUT:
%   aMat   = data matrix (npts,ntrc)
%   dt     = time sample interval [s]
%   fmin   = minimum frequency to compute FFT [Hz]
%   fmax   = maximum frequency to compute FFT [Hz]
%   df     = frequency sample interval [Hz]
% OUTPUT:
%   aMatF  = complex FFT matrix (nfreqs,ntrc)
%   fArray = frequency vector of size nfreqs [Hz]
%
% Original by Matt Haney
% Modified by Dylan Mikesell
% Last modified 23 April 2014

%--------------------------------------------------------------------------
% Dimension for allocating
[npts,nr] = size(aMat);
%--------------------------------------------------------------------------
% Range of frequencies
%--------------------------------------------------------------------------
fmax   = ceil(fmax/df)*df;  % round fmax up to a multiple of the desired frequency spacing
fmin   = floor(fmin/df)*df; % round fmin down to a multiple of the desired frequency spacing
fArray = fmin:df:fmax;    % vector of frequencies
nf     = numel(fArray);        % number of elements in frequency vector
%--------------------------------------------------------------------------
% Pad data to length of data needed for desired frequency resolution
%--------------------------------------------------------------------------
n = round(1/(df*dt)); % must be an integer
% pad or window data to achieve length n for desired frequency resolution
if (n > npts)
    aMatPad = [aMat ; zeros(n-npts,nr)];
else
    aMatPad = aMat(1:n,:);
end
%--------------------------------------------------------------------------
% FFT each trace
%--------------------------------------------------------------------------
% cut out the frequencies desired
if (floor(n/2) == ceil(n/2)) % for even samples 
    n2 = n/2;                % start sample
else                         % for odd samples
    n2 = (n-1)/2;            % start sample
end
% Fourier transform data which has been padded to obtain the desired frequency resolution
aMatF = zeros(nf,nr); % allocate
for ii=1:nr
    dum         = fftshift(fft(aMatPad(:,ii)));
    aMatF(:,ii) = dum((n2+1+(fmin/df)):(n2+1+(fmax/df)),1); % cut out the part I want
end

return
