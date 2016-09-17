function W = Hoa2MicRadialWeightingsFromMeas(l,freq,rad,measurements)

% ATTENTION : in the case where radial weightings are computed from a set
% of measured impulse responses, it is implicitely assumed that the
% corresponding sensor has an axi-symmetric directivity. It is then much
% more efficient to work in the xOz plane. The following code requires that
% the source positions are defined by a [azm elv rad] array, where azm = 0, 
% -pi/2 <= elv <= pi/2 and rad must be the same for all the sources.
% REMARK : the results will probably not be accurate at fs/2

% Useful shortcuts
sourceSphCoord   = measurements.sourceSphCoord ;
impulseResponses = measurements.impulseResponses ;
fs               = measurements.sampFreq ;
fftNb = size(impulseResponses,1) ;
nbSou = size(impulseResponses,2) ;

% Checks if the asked frequency values are availables given the measurments
if max(freq) > max(measurements.sampFreq/2)
    error('Measurement data does not allow such high frequencies') ;
end

% Checks if the measurement data allows the hoaFmt
maxOrder = floor((nbSou-1)/2) ;
if max(l) > maxOrder
    error('There is not enough measured data to handle such high orders') ;
end

% Creates a new hoaFmt for computational needs
hoaFmtPhy = GenerateHoaFmt('res2d',0,'res3d',maxOrder) ;

% Vector of the frequency values, and corresponding wave numbers
measFreq = [.1 ; (fs/fftNb : fs/fftNb : fs/2)'] ;

% Frequency responses are computed from rearranged impulse responses
[angleMin,angleInd] = min(abs(sourceSphCoord(:,2)-pi/2)) ;
[impMax,impMaxInd] = max(abs(impulseResponses(:,angleInd))) ;
freqResponses = fft(circshift(impulseResponses, ...
    [1-impMaxInd-round(rad*fs/340) 0])) ;

% Spherical harmonic coefficients at source positions
YSou = SphericalHarmonicMatrix(hoaFmtPhy, ...
    sourceSphCoord(:,1),sourceSphCoord(:,2)) ;

% Sperical harmonic coefficients at microphone position
YMic = SphericalHarmonicMatrix(hoaFmtPhy,0,pi/2) ;

% Microphone radial weightings are computed for each measurement freq.
W = zeros(hoaFmtPhy.nbComp,length(measFreq)) ;
invYsou = pinv(YSou,1e-6) ;
for I = 1 : length(measFreq)
    W(:,I) = ( ( freqResponses(I,:) * invYsou ) ./ YMic.' ).' ;
end

% Only the requested weightings are kept
W = W((0:maxOrder)<=max(l),:) ;

% Weightings are interpolated to match the desired frequency values
W = interp1(measFreq(:),W.',freq(:),'spline').' ;

% Formatting W
W = W(l+1,:) ;


