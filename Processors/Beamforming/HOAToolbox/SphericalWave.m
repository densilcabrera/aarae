function P = SphericalWave(xyzSou,xyzPts,freq,origGain)

% Check that the points and speakers are coherent in terms of dimension
if size(xyzPts,2) ~= size(xyzSou,2)
    error('Source and point position arrays must have the same number of columns') ;
end

% Default (complex) gain at the origin
if nargin < 4
    origGain = [] ;
end

% Frequency etc
freq = freq(:) ;
nbFreq = length(freq) ;
nbPts  = size(xyzPts,1) ;
nbSou  = size(xyzSou,1) ;

% Source to microphone distance matrix
Dst = XyzDistanceMatrix(xyzSou,xyzPts)' ;
Dst(Dst==0) = eps ;

% Pressure array
P = zeros(nbPts,nbSou,nbFreq) ;
if isempty(origGain)
    for I = 1 : length(freq)
        P(:,:,I) = exp( Dst * (-1i*2*pi*freq(I)/340) ) ./ Dst ;
    end
else
    Rad = repmat(sqrt(sum(abs(xyzSou).^2,2)).',nbPts,1) ;
    for I = 1 : length(freq)
        P(:,:,I) = origGain ...
             * exp( Rad * ( 1i*2*pi*freq(I)/340) ) .* Rad ...
            .* exp( Dst * (-1i*2*pi*freq(I)/340) ) ./ Dst ;
    end
end