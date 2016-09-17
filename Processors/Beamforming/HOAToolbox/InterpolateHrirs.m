function HriInt = InterpolateHrirs(azm,elv,Hri,azmInt,elvInt,smpFrq)


% Check that the hrirs are formatted appropriately
if size(Hri,2)~=2
    error(['HRIRs must be passed in an array of dimension' ...
    '[ nmbTap x 2 x nmbDir ] '])
end 
if ~ ( (size(Hri,3)==length(azm)) && (size(Hri,3)==length(elv)) )
    error(['The HRIR array is not consistent with the azimuth and ' ...
        'elevation vectors']) ;
end

% Check that azmInt is consistent with elvInt
if length(azmInt) ~= length(elvInt)
    error(['The interpolation azimuth and elevation arrays must' ...
        ' have the same number of elements']) ;
end


% default sampling frequency is 48kHz
if nargin < 6
    smpFrq = 48000 ;
end


% Length of the filters
nmbTap = size(Hri,1) ;

% number of input (measurement) directions
nmbMea = length(azm) ;

% Interpolate the HRTF magnitude using splines
Hrt = permute(abs(fft(Hri)),[1 3 2]) ;
HrtInt(:,:,1) = SphericalSplineInterp(azm,elv, ...
    Hrt(1:end/2+1,:,1).',azmInt,elvInt) ;
HrtInt(:,:,2) = SphericalSplineInterp(azm,elv, ...
    Hrt(1:end/2+1,:,2).',azmInt,elvInt) ;
HrtInt = permute(HrtInt,[2 3 1]) ;

% "No-phase" HRIRs
HrtInt(nmbTap/2+2:nmbTap,:,:) = HrtInt(nmbTap/2:-1:2,:,:) ;
HriInt = fftshift(real(ifft(HrtInt)),1) ;

% Corresponding minimum-phase filters
for I = 1 : 2
    for J = 1 : length(azmInt)
        [~,HriInt(:,I,J)] = rceps(HriInt(:,I,J)) ;
    end
end

% Calculate the delays between every original HRIR and its minimum-phase
% counterpart
del = zeros(nmbMea,2) ;
for I = 1 : nmbMea
    [rcp,HriRef] = rceps(Hri(:,1,I)) ;
    cor = xcorr(resample(Hri(:,1,I),2,1),resample(HriRef,2,1)) ;
    [maxCor,maxIdx] = max(cor) ;
    del(I,1) = (maxIdx-2*nmbTap)/2/smpFrq ;
    [rcp,HriRef] = rceps(Hri(:,2,I)) ;
    cor = xcorr(resample(Hri(:,2,I),2,1),resample(HriRef,2,1)) ;
    [maxCor,maxIdx] = max(cor) ;
    del(I,2) = (maxIdx-2*nmbTap)/2/smpFrq ;
end

% Interpolate the delays using spherical splines
delInt(:,1) = SphericalSplineInterp(azm,elv,del(:,1),azmInt,elvInt,1,1e-5) ;
delInt(:,2) = SphericalSplineInterp(azm,elv,del(:,2),azmInt,elvInt,1,1e-5) ;
delInt = delInt - min(min(delInt)) ;

% Delay the minimum-phase HRIRs according to these delay values
frq = (0:nmbTap/2)'*smpFrq/nmbTap ;
HrtInt = permute(fft(HriInt),[1 3 2]) ;
HrtInt(1:nmbTap/2+1,:,1) = HrtInt(1:nmbTap/2+1,:,1) ...
    .* exp(-1i*2*pi*frq*delInt(:,1).') ;
HrtInt(1:nmbTap/2+1,:,2) = HrtInt(1:nmbTap/2+1,:,2) ...
    .* exp(-1i*2*pi*frq*delInt(:,2).') ;
HrtInt(nmbTap/2+2:nmbTap,:,:) = conj(HrtInt(nmbTap/2:-1:2,:,:)) ;
HriInt = permute(real(ifft(HrtInt)),[1 3 2]) ;



