function Y = SphericalHarmonicMatrix(hoaFmt,azm,elv)
%
% ATTENTION : subsiste un petit probleme de normalisation dans le cas des
% harmoniques complexes (elles devraient toutes etre divisees par racine
% de 2 sauf la premiere). A terme, convCoef devrait etre calcule a l'aide
% de hoaConventionConversionMachin qui prendra en compte la chose.
% MERCI DE VOTRE COMPREHENSION - G. ABITBOL
%
% Hahaha! you can't read French!
%
% Y = ComputeSphericalHarmonics(hoaFmt,azm,elv)
%
% compute spherical harmonics values for given angular positions
%
% Input variables :
%    - hoaFmt : HOA format structure describing a set of shperical
%               harmonics. Type 'help GenerateHoaFmt' for details
%    - azm   : vector or matrix of azimuth values
%    - elv   : vector or matrix of elevation values
%               elv array must be the same size as azm array
%               elv array is optional, the defaut value is zero 
%
% Output :
%    Array of spherical harmonics values 
%    with dimension [hoaFmt.nbComp x size(azm,1) x size(azm,2)]
%
%
% N. Epain on a J. Daniel idea - Last update: 18/12/2007


% PRELIMINARY TESTS AND DEFAUT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    elv = zeros(size(azm));
elseif ~isequal(size(azm), size(elv))
    error('azimuth and elevation arrays must be the same size');
end

% COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

azmSze = size(azm) ;
azm = azm(:)' ;
elv = elv(:)' ;

% Normalization coefficients depending on HOA convention
switch hoaFmt.conv
    case 'SN3D'
        convCoef = @(l) 1 ;
    case 'N3D'
        convCoef = @(l) sqrt(2*l+1) ;
    case 'SN2D'
        convCoef = @(l) exp( 1/2*log(2*l+1) + l*log(2) ...
            + gammaln(l+1) - 1/2*gammaln(2*l+1+1) -1/2*log(2) ) ;
    case 'N2D'
        convCoef = @(l) exp( 1/2*log(2*l+1) + l*log(2) ...
            + gammaln(l+1) - 1/2*gammaln(2*l+1+1) ) ;
end

% Computation of spherical harmonics values (real or complex ones)
% The different spherical harmonics are sorted to make computations faster
Y = zeros(hoaFmt.nbComp,length(azm)) ;
[idx,prm] = sortrows(hoaFmt.index) ;
l = idx(:,1) ;
m = idx(:,2) ;
switch hoaFmt.type
    case 'real'
        J = -1 ;
        for I = 1 : hoaFmt.nbComp
            if l(I) ~= J
                Pm = convCoef(l(I)) * legendre(l(I), sin(elv), 'sch') ;
            end
            Y(I,:) = Pm(abs(m(I))+1,:) .* ...
                ( (m(I)<0 ) .* sin(-m(I)*azm) ...
                + (m(I)==0) ...
                + (m(I)>0 ) .* cos( m(I)*azm) ) ;
            J = l(I) ;
        end
    case 'complex'
        J = -1 ;
        for I = 1 : hoaFmt.nbComp
            if l(I) ~= J
                Pm = convCoef(l(I)) * legendre(l(I), sin(elv), 'sch') ;
            end
            Y(I,:) = Pm(abs(m(I))+1,:) .* ...
                ( (m(I)==0) + (m(I)~=0)/sqrt(2)*exp(1i*m(I)*azm) ) ;
            J = l(I) ;
        end
end

% Matrix Y is sorted back
Y(prm,:) = Y ;

% Matrix Y is given the desired format
Y = squeeze(reshape(Y, [size(Y,1) azmSze])) ;
