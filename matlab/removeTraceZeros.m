A = W(3);
plot(A)



clc

Ad = double(A);

Ad = [zeros(10000,1); Ad(1:100000)];
% Ad = [zeros(100000,1); Ad; Ad; Ad(1:10000); Ad(end-10000:end)];
npts = numel(Ad);

subplot(2,2,1);
plot(Ad);
title('Origianl Trace');
axis('tight');

ZeroIdx = (Ad == 0);
subplot(2,2,2);
plot(ZeroIdx);
title('Zero samples');
axis('tight');

AdZ = Ad;
AdZ(ZeroIdx) = [];
subplot(2,2,3);
plot(AdZ);
title('Trace after removing zeros');
axis('tight');
%---------------------------------------------------------------------
% Use derivative to find edges of zero blocks
dZeroIdx = diff(ZeroIdx);
startIdx = find(dZeroIdx == 1); % beginning of zero blocks
stopIdx  = find(dZeroIdx == -1); % end of zero blocks
AdR      = zeros( npts, 1 ); % allocate new trace to reconstruct

startZeroFlag = 0;
if isempty(startIdx) && isempty(stopIdx) == 0
    startIdx = 1;
    startZeroFlag = 1;
end
if isempty(stopIdx) && isempty(startIdx) == 0
    stopIdx = numel(Ad);
end
%---------------------------------------------------------------------
% Treat the FRONT of the trace
if stopIdx(1) < startIdx(1) % trace starts with zeros
    startIdx = [1; startIdx];
else
    AdR( 1 : startIdx(1) ) = AdZ( 1 : startIdx(1) );
end
%---------------------------------------------------------------------
% Treat the BACK of the trace
if stopIdx(end) < startIdx(end) % trace ends with zeros
    stopIdx = [stopIdx; npts];
else % trace has data at end
    nSamples = npts - stopIdx(end) - 1;
    AdR(end-nSamples:end) = AdZ(end-nSamples:end);
end
%---------------------------------------------------------------------
% Treat the center parts of the trace
nWindows = numel(startIdx); % number of zero windows
start    = 1; % the local starting counter for the data with no zeros

for ii = 1 : nWindows-1
    nSamples = startIdx(ii+1) - stopIdx(ii) - 1; % number of data points in this window
    AdR( stopIdx(ii) + 1 : stopIdx(ii) + 1 + nSamples ) = AdZ( start : start + nSamples );
    start = start + nSamples + 1; % update where we start next
end

if startZeroFlag
    AdR(1) = 0; % have to treat frist point
end

subplot(2,2,4);
plot(AdR);
title('Reconstructed trace with zeros');
axis('tight')

fprintf('Sum of trace difference = %0.30f\n',sum(Ad - AdR));

figure;
plot(AdR-Ad);











