function W = Hoa2MicRadialWeightings(l,kr,varargin)
% 
% Computes the "radial weightings" or "modal coefficients" Wn(kr) of the
% Fourier-Bessel expansion for various types of sensors, including sphere
% and cylinder-baffled microphones.
%
% NOTE: in the case where the baffle is a cylinder, the expansion is
% implicitely chosen as a CYLINDRICAL harmonic expansion, which means the
% spherical Bessel and Hankel functions become cylindrical Bessel and
% Hankel functions.

% Default values for optional arguments
micType = 'omni' ;
kR = zeros(size(kr)) ;
admittance = zeros(size(kr)) ;
freqResp = [] ;
baffleType = 'sphere' ;

% Number of arguments is checked
if round(length(varargin)/2)~=length(varargin)/2
    error('illegal number of arguments') ;
end

% Test and assignment of the optional arguments
for I = 1 : 2 : length(varargin)-1
    switch varargin{I}
        
        case 'micType'
            switch varargin{I+1}
                case {'omni','cardio','super','hyper','eight'}
                    micType = varargin{I+1} ;
                otherwise
                    error([varargin{I+1} ...
                         ' is not a valid microphone type']) ;
            end

        case 'kR'
            if strcmp(class(varargin{I+1}),class(kr))
                if size(varargin{I+1})==size(kr)
                    kR = varargin{I+1} ;
                    if ~isempty(kr(kr<kR))
                        if isempty(kr(kr<(1-1e-6)*kR))
                            kr(kr<(1-1e-6)*kR) = kR ;
                        else
                            error('kR cannot be larger than kr.') ;
                        end
                    end
                elseif isempty(varargin{I+1})
                else
                    error('kR must be the same size as kr or empty') ;
                end
            else
                error('kR and kr must be the same data type') ;
            end

        case 'impedance'
            if strcmp(class(varargin{I+1}),class(kr))
                if isscalar(varargin{I+1})
                    admittance = 1 / varargin{I+1} * ones(size(kr)) ;
                elseif size(varargin{I+1})==size(kr)
                    admittance = 1 ./ varargin{I+1} ;
                elseif isempty(varargin{I+1})
                else
                    error(['Impedance must be either the same size as' ...
                           ' kr or a scalar, or empty ']);
                end
            else
                error('Impedance must be the same data type as kr') ;
            end
            
        case 'freqResp'
            if strcmp(class(varargin{I+1}),class(kr))
                if size(varargin{I+1}) == size(kr)
                    freqResp = varargin{I+1} ;
                elseif isempty(varargin{I+1}) ; 
                else
                    error(['Frequency response must be either the same' ...
                        ' size as kr or empty ']);
                end
            else
                error('Frequency response must be the same data type as kr') ;
            end

        case 'baffleType'
            switch varargin{I+1}
                case {'sphere','cylinder'}
                    baffleType = varargin{I+1} ;
                otherwise
                    error([varargin{I+1} ...
                         ' is not a valid baffle type']) ;
            end

        otherwise
            error(['''' varargin{I} '''... what''s that??']) ;
    end
end

% kr, l, kR and admittance are formatted
krSize = size(kr) ;
kr = kr(:)' ;
[kr,l] = meshgrid(kr,l) ;
kR = repmat(kR(:)',size(l,1),1) ;
admittance = repmat(admittance(:)',size(l,1),1) ;

% Coefficients a & b (see below) depending on the type of sensors
switch micType
    case 'omni'
        a = 1 ; 
    case 'cardio'
        a = 1/2 ;
    case 'super'
        a = 3/8 ;
    case 'hyper'
        a = 1/4 ;
    case 'eight'
        a = 0 ;
end
b = 1 - a ;

% Weightings corresponding to the direct sound field
switch baffleType
    
    case 'sphere'
        
        W = 1i.^l ...
            .* (a*SphericalBesselj(l,kr)-1i*b*SphericalBesseljDerivative(l,kr)) ;
        
    case 'cylinder'
        
        W = 1i.^l ...
            .* (a*besselj(l,kr)-1i*b*BesseljDerivative(l,kr)) ;
        
end

% Weightings corresponding to the diffracted soundfield (if any)
if any(any(kR))
    
    switch baffleType
    
        case 'sphere'

            W(kr~=0) = W(kr~=0) - 1i.^l(kr~=0) ...
                .* (a*SphericalBesselh(l(kr~=0),2,kr(kr~=0))-1i*b*SphericalBesselhDerivative(l(kr~=0),2,kr(kr~=0))) ...
                .* (SphericalBesseljDerivative(l(kr~=0),kR(kr~=0))+1i*admittance(kr~=0).*SphericalBesselj(l(kr~=0),kR(kr~=0))) ...
                ./ (SphericalBesselhDerivative(l(kr~=0),2,kR(kr~=0))+1i*admittance(kr~=0).*SphericalBesselh(l(kr~=0),2,kR(kr~=0))) ;
            
        case 'cylinder'

            W(kr~=0) = W(kr~=0) - 1i.^l(kr~=0) ...
                .* (a*besselh(l(kr~=0),2,kr(kr~=0))-1i*b*BesselhDerivative(l(kr~=0),2,kr(kr~=0))) ...
                .* (BesseljDerivative(l(kr~=0),kR(kr~=0))+1i*admittance(kr~=0).*besselj(l(kr~=0),kR(kr~=0))) ...
                ./ (BesselhDerivative(l(kr~=0),2,kR(kr~=0))+1i*admittance(kr~=0).*besselh(l(kr~=0),2,kR(kr~=0))) ;
            
    end

end

% Apply the frequency response if necessary
if ~isempty(freqResp)
    W = W .* repmat(freqResp(:).',size(l,1),1) ;
end

% Matrix W is given the desired format
W = squeeze(reshape(W, [size(W,1) krSize])) ;

% Some values can be NaN for l very large and kr very small.
% These NaNs are replaced by zeros
W(isnan(W)) = 0 ;

end