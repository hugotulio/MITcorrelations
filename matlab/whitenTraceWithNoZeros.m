function Aout = whitenTraceWithNoZeros(A, F, dt)
%
% This function spectrally whitens the trace A.
%
% USAGE: Aout = whitenTrace(A, F, dt)
%
% INPUT:
%   A  = the trace (npts,1)
%   F  = a vector of min and max frequencies [fmin,fmax]
%   dt = sample interval in seconds
% OUTPUT:
%   Aout = the whitened trace (npts,1)
% NOTE:
%   When reinserting the zero blocks, no tapering is applied. This
%   does not appear to cause any problems, but would be useful to
%   taper 5% of each edge or something like this. This would only be a
%   minor improvement.
%
% Written by Piero Poli (Original was blanchmat.m)
% Modified by Dylan Mikesell (mikesell@mit.edu)
% Last modified 2 June 2014

npts = numel(A);
Aout = zeros( npts, 1 ); % allocate new trace to reconstruct

% check for zeros
zeroIdx = (A == 0); % indices of zero traces

%---------------------------------------------------------------------
% Remove parts of trace that are zero
if sum(zeroIdx) > 0
        
    A(zeroIdx) = []; % remove zeros
    Awhite     = whitenTrace(A, F, dt); % whiten the data

    %---------------------------------------------------------------------
    % Instert parts of trace that are zero
    
    %---------------------------------------------------------------------
    % Use derivative to find edges of zero blocks
    dZeroIdx = diff( zeroIdx ); % derivative of
    startIdx = find( dZeroIdx ==  1 ); % beginning of zero blocks
    stopIdx  = find( dZeroIdx == -1 ); % end of zero blocks

    %-----------------------------------------------------------------
    % Treat FRONT and BACK if only 1 zero block
    startZeroFlag = 0;
    if isempty(startIdx) && isempty(stopIdx) == 0 % trace starts with zeros with only one zero block.
        startIdx = 1;
        startZeroFlag = 1;
    end
    if isempty(stopIdx) && isempty(startIdx) == 0 % trace ends with zeros with only one zero block.
        stopIdx = npts;
    end
    %---------------------------------------------------------------------
    % Treat the FRONT of the trace
    if stopIdx(1) < startIdx(1) % trace starts with zeros
        startIdx = [1; startIdx];
    else
        Aout( 1 : startIdx(1) ) = Awhite( 1 : startIdx(1) );
    end
    %---------------------------------------------------------------------
    % Treat the BACK of the trace
    if stopIdx(end) < startIdx(end) % trace ends with zeros
        stopIdx = [stopIdx; npts];
    else % trace has data at end
        nSamples = npts - stopIdx(end) - 1;
        Aout(end-nSamples:end) = Awhite(end-nSamples:end);
    end
    %---------------------------------------------------------------------
    % Treat the center parts of the trace
    nWindows = numel(startIdx); % number of zero windows
    start    = 1; % the local starting counter for the data with no zeros
    
    for ii = 1 : nWindows-1
        nSamples = startIdx(ii+1) - stopIdx(ii) - 1; % number of data points in this window
        Aout( stopIdx(ii) + 1 : stopIdx(ii) + 1 + nSamples ) = Awhite( start : start + nSamples );
        start = start + nSamples + 1; % update where we start next
    end
    
    if startZeroFlag
        Aout(1) = 0; % have to treat frist point
    end
    
else % no zero blocks so whiten as normal
    
    Aout = whitenTrace(A, F, dt); % whiten the data
    
end

return