function Trn = PlaneWave(xyzSou,xyzPts,frq)

% Test the dimension of the spaces where the sources and points are defined
dimSou = size(xyzSou,2) ;
dimPts = size(xyzPts,2) ;
if dimSou ~= dimPts
    error('The source and point position arrays must have the same number of columns') ;
end

% Shortcuts
dim    = dimSou ;
frq    = frq(:) ;
nmbFrq = length(frq) ;
nmbPts = size(xyzPts,1) ;
nmbSou = size(xyzSou,1) ;

% Unitary vectors pointing towards the sources
dirSou = -xyzSou ./ repmat(sqrt(sum((xyzSou).^2,2)),1,dim) ; 

% Dot-product of the point vectors with the source direction vectors
dotPrd = xyzPts * dirSou';

% Array of the frequency transfers between the sources and the points
Trn = zeros(nmbPts,nmbSou,nmbFrq) ;
for I = 1 : length(frq)
    Trn(:,:,I) = exp( dotPrd * (-1i*2*pi*frq(I)/340) ) ;
end

end
