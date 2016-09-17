function order = Hoa2MicRequiredOrder(micFmt,freq,minOrder,maxOrder,tolerance)

% Default tolerance
if nargin < 5
    tolerance = 1e-3 ;
end

% Default maximum order
if nargin < 4
    maxOrder = 100 ;
end

% Default minimum order
if nargin < 3
    minOrder = max(floor(sqrt(micFmt.nbMic)-1),1) ;
end


% Wavenumber
k = 2*pi*max(freq(:))/340 ;

% The order required for an accurate modelling is computed
minOrderReq = minOrder ;
maxOrderPos = maxOrder + 1 ;
for I = 1 : micFmt.nbArrays
    switch micFmt.arrays(I).micType
        case 'measured'
            nbPos = size(micFmt.arrays(I).measurements.sourceSphCoord,1) ; 
            if size(micFmt.arrays(I).measurements.impulseResponses,3) == 1
                maxOrderPos = min(floor(nbPos/2-1),maxOrderPos) ;
                W = Hoa2MicRadialWeightingsFromMeas(0:maxOrderPos, ...
                    max(freq(:)),micFmt.arrays(1).sphCoord(1,3), ...
                    micFmt.arrays(1).measurements) ;
                minOrderReq = ...
                    max(length(find(abs(W)/max(abs(W))>tolerance))-1, ...
                    minOrderReq) ;
            else
                minOrderReq = max(round(sqrt(nbPos)-2),minOrderReq) ;
                maxOrderPos = min(round(sqrt(nbPos)-2),maxOrderPos) ;
            end
        otherwise
            W = Hoa2MicRadialWeightings(0:maxOrder, ...
                k*max(micFmt.arrays(I).sphCoord(:,3)), ...
                'micType',micFmt.arrays(I).micType,...
                'kR',k*micFmt.sphereRadius, ...
                'impedance',micFmt.sphereImpedance(freq)) ;
            minOrderReq = ...
                max(length(find(abs(W)/max(abs(W))>tolerance))-1, ...
                minOrderReq) ;
    end
end
order = min(minOrderReq,maxOrderPos) ;

if order < minOrderReq
    fprintf([ '\n WARNING: '...
        'The maximum spherical harmonic order (' num2str(order) ...
        ') has been imposed by the measurement data \n']) ;
elseif order == maxOrder
    fprintf([ '\n WARNING: '...
        'Maximum spherical harmonic order (' num2str(order) ...
        ') has been reached, consider using a higher maxOrder \n']) ; 
end
