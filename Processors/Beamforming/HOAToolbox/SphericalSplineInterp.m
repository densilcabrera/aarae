function DatOut = SphericalSplineInterp(azmInp,elvInp,DatInp,azmOut,elvOut,ord,smo)
%
%   DatOut = SphericalSplineInterp(azmInp,elvInp,DatInp,azmOut,elvOut)
%
% Interpolate input values to the desired output directions on the sphere 
% using Spherical Thin Plate Splines (STPS).
% See : G. Wahba. " Spline interpolation and smoothing on the sphere". 
%       (SIAM J. Sci. Stat. Comp., 2:5–16, 1981).
% 
% OUTPUT:
%     - DatOut: Interpolated data, at the directions specified 
%               in azmOut and elvOut.
%
% INPUTS:
%     - azmInp, elvInp:  input azimuth and elevation vectors (in rad)
%     - DatInp: Input values at directions specified in azmInp and elvInp.
%     - azmOut, elvOut: output azimuth and elevation vectors (in rad)
%
% OPTIONAL PARAMETERS: 
%     - ord: order of the STPS interpolation (see Wahba, 1981).
%            Default value is 1.
%     - smo: smoothing parameter (see Wahba, 1981).
%            lambda = 0 : no smoothing, pure interpolation
%            The effect of lambda is increasing ~exponentially. 
%            (for HRTF log-magnitude, trry with O(1e-5) first).
%            Default value is 1e-6.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pierre GUILLON, Carlab 2010, 
%  
%   modified by N. Epain, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default smoothing parameter value
if nargin < 7
    smo = 1e-6 ;
end

% Default order value
if nargin < 6
    ord = 1 ;
end

%%% common parameters only depending on measurement directions

nmbInp = length(azmInp) ;
nmbCol = size(DatInp,2) ;

[XInp YInp ZInp] = sph2cart(azmInp,elvInp,1);
D = XInp*XInp'+YInp*YInp'+ZInp*ZInp'; %dot product between input directions
D = min(D,1);
W = (1-D)/2;
A = zeros(size(D));
A(D~=1) = log(1+1./sqrt(W(D~=1)));
C = 2*sqrt(W);

if ord == 1
    Q=2*A.*W-C+1;
    Rn=(Q-.5)/(2*pi);
elseif ord == 2
    Q= (A.*(12*W.^2-4*W)-6*C.*W+6*W+1)/2;
    Rn= (.5*Q-1/6)/(2*pi);  
else
    error('Error : ord should be equal to 1 or 2');
end
T = ones(nmbInp,1);
M = Rn+nmbInp*smo*eye(nmbInp);
Minv = M^(-1) ;
    
%%% parameters depending on measured data
c = zeros(nmbInp,nmbCol) ;
d = zeros(nmbCol,1) ;
for I = 1:nmbCol    
    d(I,1) = (-T'*Minv*DatInp(:,I))/(T'*Minv*T);
    c(:,I) = Minv*(DatInp(:,I)+d(I,1)*T);
end

%%% interpolation to target directions
nmbOut = length(azmOut) ;
DatOut = zeros(nmbOut,nmbCol) ;
[XOut YOut ZOut] = sph2cart(azmOut,elvOut,1);
for J = 1:nmbOut

    D = XInp*XOut(J)+YInp*YOut(J)+ZInp*ZOut(J);
    D = min(D,1);
    W = (1-D)/2;
    A = zeros(size(D));
    A(D~=1) = log(1+1./sqrt(W(D~=1)));
    C = 2*sqrt(W);

    if ord==1

        Q = 2*A.*W-C+1;
        Rt = (Q-.5)/(2*pi);

    elseif ord==2

        Q = (A.*(12*W.^2-4*W)-6*C.*W+6*W+1)/2;
        Rt = (.5*Q-1/6)/(2*pi);

    end

    DatOut(J,:) = Rt.'*c-d.' ;

end