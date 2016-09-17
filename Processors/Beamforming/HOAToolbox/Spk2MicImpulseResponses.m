function impRsp = Spk2MicImpulseResponses(spkFmt,micFmt,impLng,smpFrq)

% Default sampling frequency
if nargin < 4
    smpFrq = 48000 ;
end

% Default length for the impulse responses
if nargin < 3
    impLng = 256 ;
end

% Are the impulse responses long enough ?
xyzMic = zeros(micFmt.nbMic,3) ;
sphMic = zeros(micFmt.nbMic,3) ;
micCnt = 0 ;
for I = 1 : micFmt.nbArrays
    xyzMic(micCnt+1:micCnt+micFmt.arrays(I).nbMic,:) = ...
        micFmt.arrays(I).xyzCoord ;
    sphMic(micCnt+1:micCnt+micFmt.arrays(I).nbMic,:) = ...
        micFmt.arrays(I).sphCoord ;
    micCnt = micCnt + micFmt.arrays(I).nbMic ;
end
if ~any(isinf(spkFmt.sphCoord(:,3)))
    % Then there are only spherical sources.
    % What is the difference between the max and min spk to mic delay ?
    xyzSpk = spkFmt.xyzCoord ;
    SpkMicDst = XyzDistanceMatrix(xyzSpk,xyzMic) ;
    maxDel = round((max(max(SpkMicDst))-min(min(SpkMicDst)))/340*smpFrq) ;
elseif ~any(~isinf(spkFmt.sphCoord(:,3)))
    % Then there are only plane-wave sources.
    % What is the maximum mic to mic delay ?
    avgMic = mean(xyzMic,1) ;
    [maxDst,farIdx] = max(sphMic(:,3)) ;
    farMic = xyzMic(farIdx,:) ;
    micMicDst = XyzDistanceMatrix(farMic,avgMic) ;
    maxDel = round(max(max(2*micMicDst))/340*smpFrq) ;
else
    % It won't work
    error(['The spkFmt contains plane-waves as well as spherical ' ...
        'sources. This is forbidden as it leads to infinite delays']) ;
end
if impLng < 2*maxDel
     warning('HOA:ImpulseResponsesAreTooShort', ...
         'Consider using larger impulse responses') ;
end

% Number of sources and microphones
nmbSpk = spkFmt.nbSpk ;
nmbMic = micFmt.nbMic ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPULSE RESPONSE CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FFT number of points
fftNmb = 2*impLng ;

% FFT frequency values and corresponding wave numbers
fftFrq = [min(1,smpFrq/fftNmb/2),smpFrq/fftNmb:smpFrq/fftNmb:smpFrq/2]' ;

% Transfer matrix at every frequency. The gains and delays are chosen so
% that the average Spk -> Mic distance corresponds to a unit gain with no
% delay. This is done by using the 'reference radius' option is the
% transfer routine.
if ~any(isinf(spkFmt.sphCoord(:,3)))
    avgMic = mean(xyzMic,1) ;
    SpkAvgMicDst = XyzDistanceMatrix(xyzSpk,avgMic) ;
    radRef = mean(SpkAvgMicDst) ;
else
    radRef = [] ;
end
impRsp = Spk2MicTransfers(spkFmt,micFmt,fftFrq,radRef) ;
impRsp = permute(impRsp,[3 1 2]) ;

% Back in the time domain
impRsp(fftNmb/2+2:fftNmb,:,:) = conj(impRsp(fftNmb/2:-1:2,:,:)) ;
impRsp = fftshift(real(ifft(impRsp)),1) ;

% Time-windowing
wdw = hanning(impLng/2) ;
wdw = [wdw(1:end/2);ones(impLng/2,1);wdw(end/2+1:end)] ;
impRsp = repmat(wdw,[1 nmbMic nmbSpk]) ...
    .* impRsp(1/2*impLng+1:3/2*impLng,:,:) ;
