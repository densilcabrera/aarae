function W = Spk2HoaRadialWeightings(l,kr)

% kr and l are formatted
krSize = size(kr) ;
kr = kr(:)' ;
plwIdx = isinf(kr) ;
[kr,l] = meshgrid(kr,l) ;

% Weightings
W = 1i.^(-l) .* SphericalBesselh(l,2,kr) ./ SphericalBesselh(0,2,kr) ;

% Plane-wave case (infinite kr)
W(:,plwIdx) = 1 ;

% Matrix W is given the desired format
W = squeeze(reshape(W, [size(W,1) krSize]));

end