function D = Hoa2SpkDecodingMatrix(hoaFmt,azmSpk,elvSpk,varargin)

% default option values
invMethod = 'pinv' ;
decodType = 'basic' ;
normalise = 'none' ;

% Number of optional arguments is checked
if round(length(varargin)/2)~=length(varargin)/2
    error('illegal number of arguments') ;
end

% New option values are checked and assigned
for I = 1 : 2 : length(varargin)-1
    switch varargin{I}
        case 'invMethod'
            switch varargin{I+1}
                case {'transpose','pinv','tikhonov'}
                    invMethod = varargin{I+1} ;
                otherwise
                    error(['Available inversion methods are ''pinv'' ' ...
                        'and ''transpose''']) ;
            end
        case 'decodType'
            switch varargin{I+1}
                case {'basic','maxrE','inPhase'}
                    decodType = varargin{I+1} ;
                otherwise
                    error(['Available decoding schemes are ''basic'' ' ...
                        ', ''maxrE'', or ''inPhase''']) ;
            end
        case 'normalise'
            switch varargin{I+1}
                case {'none','energy'}
                    normalise = varargin{I+1} ;
                otherwise
                    error(['Available normalisations are ''none'' ' ...
                        'and ''energy''']) ;
            end
        otherwise
            error('Unknown option') ;
    end
end

% Shortcuts
nmbHrm = hoaFmt.nbComp ;
nmbSpk = length(azmSpk) ;
maxOrd = max(hoaFmt.index(:,1)) ;
switch hoaFmt.conv
    case {'N2D','SN2D'}
        decDim = '2d' ;
    case {'N3D','SN3D'}
        decDim = '3d' ;
end

% Compute the matrix of spherical harmonic coefficients
Y = SphericalHarmonicMatrix(hoaFmt,azmSpk,elvSpk) ;

% Compute the gain matrix
switch invMethod
    case 'pinv'
        D = pinv(Y) ;
    case 'tikhonov'
        if nmbSpk >= nmbHrm
            D = Y'/(Y*Y'+1e-2*eye(nmbHrm)) ;
        else
            D = (Y'*Y+1e-2*eye(nmbSpk))\Y' ;
        end
    case 'transpose'
        A = zeros(nmbHrm^2,nmbSpk) ;
        d = zeros(nmbHrm^2,1) ;
        for I = 1 : nmbHrm
            for J = 1 : nmbHrm
                A((I-1)*nmbHrm+J,:) = Y(I,:).*Y(J,:) ;
                d((I-1)*nmbHrm+J) = (I==J) ;
            end
        end
        wgt = pinv(A) * d ;
        D = diag(wgt) * Y.' ;
 end

% Order weightings
switch decodType
    case 'maxrE'
        weights = Hoa2SpkMaxrEWeightings(maxOrd,decDim) ;
        Wgt = sparse(1:nmbHrm,1:nmbHrm,weights(hoaFmt.index(:,1)+1)) ;
        D = D * Wgt ;
        switch normalise
            case 'energy'
                D = D * sqrt(nmbHrm) / norm(Wgt,'fro') ;
        end
    case 'inPhase'
        weights = Hoa2SpkInPhaseWeightings(maxOrd,decDim) ;
        Wgt = sparse(1:nmbHrm,1:nmbHrm,weights(hoaFmt.index(:,1)+1)) ;
        D = D * Wgt ;
        switch normalise
            case 'energy'
                D = D * sqrt(nmbHrm) / norm(Wgt,'fro') ;
        end
end

