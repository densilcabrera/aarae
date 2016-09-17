function G = Vbap3dGainMatrix(sphSpk,sphSou,vectOpt,normOpt,spkTri)
%
%  Calculate a VBAP panning gain matrix for a 3D loudspeaker setup 
%
%  G = Vbap3dGainMatrix(sphSpk,sphSou,[vectOpt],[normOpt],[spkTri])
%
%  Inputs:
%  - sphSpk contains the azimuths and elevations of the speakers [Mx2]
%  - sphSou contains the azimuths and elevations of the sources [Nx2]
%  - vectOpt (optional) sets the "vector" mode of the panning. When set to
%    'velocity' (default, classic VBAP method) the panning gains are 
%    calculated so that the velocity vector is in the source direction. 
%    When set to 'energy', the gains ensure the energy vector is in the
%    source direction.
%  - spkTri (optional) contains the triangulation (convex hull) of the
%    speaker points. This saves some time when calculating some panning 
%    gains for different source positions with the same speaker setup in
%    a loop. If set to empty (default) the convex hull is calculated by
%    the routine.
%  - normOpt (optional) determines the way the gains are normalised.
%    If set to 'amplitude' (default), the sum of the directional gains is 1 
%    for every direction. If set to 'energy' then the sum of the squared
%    gains is 1 for every source direction. 
%
%  Output:
%  - G is a real matrix of dimension [MxN]


% Default value for the speaker triangulation
% In the case where the function is called many times with the same speaker
% setup, some time can be saved by passing the triangulation... 
if nargin < 5 
    spkTri = [] ;
end

% Default normalization
if nargin < 4 
    normOpt = 'energy' ;
end

% Default "vector" option
if nargin < 3
    vectOpt = 'velocity' ; 
end

% Number of sources and speakers
nbSpk = size(sphSpk,1) ;
nbSou = size(sphSou,1) ;

% Cartesian coordinates of the sources and speakers (radius=1)
[xyzSpk(:,1),xyzSpk(:,2),xyzSpk(:,3)] = ...
    sph2cart(sphSpk(:,1),sphSpk(:,2),1) ;
[xyzSou(:,1),xyzSou(:,2),xyzSou(:,3)] = ...
    sph2cart(sphSou(:,1),sphSou(:,2),1) ;

% Speaker array triangulation 
if isempty(spkTri)
    spkTri = convhulln(xyzSpk) ;
end

% Finding which triangle the sources are belonging to
triInd = tsearchn([xyzSpk; 0 0 0], ...
        [spkTri (nbSpk+1)*ones(size(spkTri,1),1)],.01*xyzSou) ;

% Speaker gains
G = zeros(nbSpk,nbSou) ;
for I = 1 : nbSou
    
    % Index of the closest speakers
    spkInd = spkTri(triInd(I),:) ;
    
    % Gains for the closest speakers
    G(spkInd,I) = xyzSou(I,:) / xyzSpk(spkInd,:) ;
        
end

% Energy-vector VBAP    
if strcmpi(vectOpt,'energy')
    G = sqrt(max(G,0)) ;
end

% Normalisation: 
switch normOpt
    case 'amplitude'
        G = G * diag(1./sum(G)) ;
    case 'energy'
        G = G * diag(1./sqrt(sum(G.^2))) ;
end

