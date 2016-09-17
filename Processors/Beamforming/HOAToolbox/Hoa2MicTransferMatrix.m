function T = Hoa2MicTransferMatrix(hoaFmt,micFmt,freq)

% Maximum order of the spherical harmonics in hoaFmt
maxOrder = max(hoaFmt.index(:,1)) ;

% ShortCuts
R      = micFmt.sphereRadius ;
Z      = micFmt.sphereImpedance ;
nbComp = hoaFmt.nbComp ;
l      = hoaFmt.index(:,1) ;
m      = hoaFmt.index(:,2) ;

% Frequency vector, number of frequency values and wave number vector
freq = freq(:) ;
nbFreq = length(freq) ;
k = 2*pi*freq/340 ;

% Matrix T value is initialised
T = zeros(micFmt.nbMic,nbComp,nbFreq) ;

% Each sub-array is processed one after the other
micCnt = 0 ;
for I = 1 : micFmt.nbArrays
    
    % Shortcuts
    nbMic        = micFmt.arrays(I).nbMic ;
    micType      = micFmt.arrays(I).micType ;
    micFreqResp  = micFmt.arrays(I).micFreqResp ;
    sphCoord     = micFmt.arrays(I).sphCoord ;
    sphDirect    = micFmt.arrays(I).sphDirect ;
    measurements = micFmt.arrays(I).measurements ;
    azm          = sphCoord(:,1) ;
    elv          = sphCoord(:,2) ;
    r            = sphCoord(:,3) ;
    azmDir       = sphDirect(:,1) ;
    elvDir       = sphDirect(:,2) ;
    
    % Matrix of the spherical harmonic function values
    Y = SphericalHarmonicMatrix(hoaFmt,azm,elv).' ;
    
    % Condition : are the microphones defined by 3d antenna measurements ?
    meaCondition = false ;
    if strcmp(micType,'measured')
        meaCondition = size(measurements.impulseResponses,3) > 1 ;
    end
    
    % Condition : are the microphones are centrifugally oriented ?
    dirCondition = ~any(any(sphCoord(:,1:2)~=sphDirect)) ;
    
    % Condition : are the microphones at the same distance from the origin?
    radCondition = ~any(r~=r(1)) ;
    
    % Depending on the conditions...
    if meaCondition
        
        % Shortcuts
        sourceSphCoord   = measurements.sourceSphCoord ;
        impulseResponses = measurements.impulseResponses ;
        fs               = measurements.sampFreq ;
        nbMic = size(impulseResponses,2) ;
        nbSou = size(impulseResponses,3) ;
        
        % Checks if the asked frequency values are available given the measurments
        if max(freq) > fs/2
            error('Measurement data does not allow such high frequencies') ;
        end
        
        % Checks if the measurement data allows the chosen hoaFmt
        maxOrder = round(sqrt(nbSou)-2) ;
        if max(l) > maxOrder
            warning(['There is not enough measured data ' ...
                'to handle this hoaFmt']) ;
        end
        
        % Microphone array frequency responses
        freqResponses = zeros(nbMic,nbSou,nbFreq) ;
        for J = 1 : nbMic
            for K = 1 : nbSou
                freqResponses(J,K,:) = permute( ...
                    freqz(impulseResponses(:,J,K),1,freq,fs),[2 3 1]) ;
            end
        end
        
        % Spherical harmonic coefficients in source directions
        YSou = SphericalHarmonicMatrix(hoaFmt, ...
            sourceSphCoord(:,1),sourceSphCoord(:,2)) ;
        
        % Calculate Spherical Transfer matrix
        if all(isinf(sourceSphCoord(:,3))) ; 
            % Then the sources are plane waves
            % YSou pseudo-inverse matrix
            invYSou = YSou'/(YSou*YSou'+1e-2*eye(nbComp)) ;
            % Calculate Spherical Transfer matrix
            Ttmp = zeros(nbMic,nbComp,nbFreq) ;
            for J = 1 : nbFreq
                Ttmp(:,:,J) = freqResponses(:,:,J) * invYSou ;
            end
        else
            % There are some spherical sources
            if all(sourceSphCoord(:,3)==sourceSphCoord(1,3)) ; 
                % Then the sources are spherical and at the same distance
                % Source modal coefficients
                WSou = Spk2HoaRadialWeightings(0:max(l), ...
                    2*pi/343*sourceSphCoord(1,3)*freq) ;
                % YSou pseudo-inverse matrix
                invYSou = YSou'/(YSou*YSou'+1e-2*eye(nbComp)) ;
                % Calculate Spherical Transfer matrix
                Ttmp = zeros(nbMic,nbComp,nbFreq) ;
                for J = 1 : nbFreq
                    Ttmp(:,:,J) = freqResponses(:,:,J) * bsxfun(@times, ...
                        invYSou,1./WSou(l+1,J).') ;
                end
            else
                % Then it's a mess
                % Source modal coefficients
                WSou = Spk2HoaRadialWeightings(0:max(l), ...
                    2*pi/343*sourceSphCoord(:,3)*freq') ;
                % Calculate Spherical Transfer matrix
                Ttmp = zeros(nbMic,nbComp,nbFreq) ;
                for J = 1 : nbFreq
                    % Source Spherical Harmonic Expansion
                    BSou = YSou .* WSou(l+1,:,J) ;
                    % Source Spherical Harmonic Expansion's pseudo inverse
                    invBSou =  BSou'/(BSou*BSou'+1e-2*eye(nbComp)) ;
                    % Spherical Transfer matrix
                    Ttmp(:,:,J) = freqResponses(:,:,J) * invBSou ;
                end
            end
        end
        
        % Calculate time shift 
        [maxAmp,maxAmpIdx] = max(abs(impulseResponses)) ;
        tmeShf = (mean(mean(maxAmpIdx))-1)/fs ;
        
        % Apply time shift
        Ttmp = bsxfun(@times,Ttmp, ...
            permute(exp(1i*2*pi*freq*tmeShf),[2 3 1])) ;
                
        % Assign T in output array
        T(micCnt+1:micCnt+nbMic,:,:) = Ttmp ;
        
        % New micCnt value (current total number of microphones)
        micCnt = micCnt + nbMic ;
        
        % Go to the next sub-array
        continue
        
    end
    
    if  dirCondition
        
        if radCondition
            
            % If all the radius are identical and the microphones are
            % oriented centrifugally, then T = YW and the weights are
            % the same for every microphone
            
            % Microphone radial weightings are computed at every freq value
            if strcmp(micType,'measured')
                W = Hoa2MicRadialWeightingsFromMeas((0:maxOrder)',freq, ...
                    sphCoord(1,3),measurements) ;
            else
                W = Hoa2MicRadialWeightings((0:maxOrder)',k*r(1), ...
                    'micType',micType,'kR',k*R,'impedance',Z(freq), ...
                    'freqResp',micFreqResp(freq)) ;
            end
            if find(isnan(W))
                W(isnan(W)) = 0 ;
            end
            
            % Filling matrix T
            for J = 1 : nbFreq
                T(micCnt+1:micCnt+nbMic,:,J) = ...
                    Y * sparse(1:nbComp,1:nbComp,W(l+1,J)) ;
            end
            
        else
            
            % If the microphones are oriented centrifugally but the radius
            % are not identical then the weights have to be computed for
            % every radius value. T = YW is still true (blockwise).
            
            % Microphone radial weightings are computed at every freq value
            % for every 'r' value
            if strcmp(micType,'measured')
                W = zeros(maxOrder+1,nbMic) ;
                for J = 1 : micFMt
                    W(:,J,:) = permute(Hoa2MicRadialWeightingsFromMeas( ...
                        (0:maxOrder)',freq,sphCoord(J,3), ...
                        measurements),[1 3 2]) ;
                end
            else
                W = Hoa2MicRadialWeightings((0:maxOrder)',r*k.', ...
                    'micType',micType,'kR',R*ones(nbMic,1)*k.', ...
                    'impedance',ones(nbMic,1)*Z(freq.'), ...
                    'freqResp',ones(nbMic,1)*micFreqResp(freq.')) ;
            end
            
            % Filling matrix T
            for J = 1 : nbFreq
                T(micCnt+1:micCnt+nbMic,:,J) = Y .* W(l+1,:,J).' ;
            end
            
        end
        
        % New micCnt value (current total number of microphones)
        micCnt = micCnt + nbMic ;
        
        % Go the next sub-array
        continue
        
    end
    
    % If the microphones are not oriented centrifugally, then it is
    % much more complicated : T has to be computed element-wise.
    
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
    
    % Matrix of the kr values
    krMat = k * r.' ;
    
    % Values of the spherical bessel functions of the first kind
    % (jl(kr))
    Jkr = SphericalBesselj( ...
        repmat((0:maxOrder).',1,nbMic*nbFreq), ...
        repmat(krMat(:).',maxOrder+1,1)) ;
    Jkr = reshape(Jkr,maxOrder+1,nbFreq,nbMic) ;
    Jkr = permute(Jkr,[3 1 2]) ;
    
    % Values of the spherical bessel functions of the first kind
    % derivatives (jl'(kr))
    JkrDer = SphericalBesseljDerivative( ...
        repmat((0:maxOrder).',1,nbMic*nbFreq), ...
        repmat(krMat(:).',maxOrder+1,1)) ;
    JkrDer = reshape(JkrDer,maxOrder+1,nbFreq,nbMic) ;
    JkrDer = permute(JkrDer,[3 1 2]) ;
    
    % uRad, uAzm and uElv values
    Urad = cos(elv) .* cos(elvDir) .* cos(azm-azmDir) ...
        + sin(elv) .* sin(elvDir) ;
    Uazm = cos(elvDir) .* sin(azm-azmDir) ;
    Uelv = - sin(elv) .* cos(elvDir) .* cos(azm-azmDir) ...
        + cos(elv) .* sin(elvDir) ;
    Urad = repmat(Urad,1,nbComp) ;
    Uazm = repmat(Uazm,1,nbComp) ;
    Uelv = repmat(Uelv,1,nbComp) ;
    
    % Matrix of the Y_l^(-m) function values
    hoaFmt_ = GenerateHoaFmt('index',[l,-m]) ;
    permMat = repmat(permute(hoaFmt.index,[1 3 2]),1,nbComp) ...
        == repmat(permute(hoaFmt_.index,[3 1 2]),nbComp,1) ;
    permMat = permMat(:,:,1) & permMat(:,:,2) ;
    [perm,perm_] = find(permMat) ; %#ok<NASGU>
    Y_ = Y(:,perm) ;
    
    % Normalization coefficients depending on HOA convention
    switch hoaFmt.conv
    case 'SN3D'
        convCoef = ones(size(l)) ;
    case 'N3D'
        convCoef = sqrt(2*l+1) ;
    case 'SN2D'
        convCoef = exp( 1/2*log(2*l+1) + l*log(2) ...
            + gammaln(l+1) - 1/2*gammaln(2*l+1+1) -1/2*log(2) ) ;
    case 'N2D'
        convCoef = exp( 1/2*log(2*l+1) + l*log(2) ...
            + gammaln(l+1) - 1/2*gammaln(2*l+1+1) ) ;
    end
    convCoef = repmat(convCoef.',nbMic,1) ;
    
    % Matrix 'Truc'
    % giving the azimuthal dependency as a function of 'm'
    switch hoaFmt.type
        case 'real'
            Truc = (repmat(m.',nbMic,1)==0) ...
                +  (repmat(m.',nbMic,1)< 0) ...
                .* sin(repmat(-m.',nbMic,1).*repmat(azm,1,nbComp)) ...
                +  (repmat(m.',nbMic,1)> 0) ...
                .* cos(repmat( m.',nbMic,1).*repmat(azm,1,nbComp)) ;
            
        case {'complex','complex+'}
            Truc = exp(-1i*repmat(m.',nbMic,1).*repmat(azm,1,nbComp)) ;
    end
    
    % Values of the associated legendre function derivatives
    [index,perm] = sortrows(hoaFmt.index) ;
    Pder = zeros(nbMic,nbComp) ;
    deg = 0 ;
    PlDer = AssociatedLegendreDerivative(0,sin(elv)) ;
    for J = 1 : nbComp
        if index(J,1) ~= deg
            PlDer = AssociatedLegendreDerivative(index(J,1),sin(elv)) ;
        end
        Pder(:,J) = PlDer(abs(index(J,2))+1,:).' ;
        deg = index(J,1) ;
    end
    Pder(:,perm) = Pder ;
    
    % A block of the Matrix T is computed at every frequency value
    for J = 1 : nbFreq
        
        % Part of the matrix corresponding to the direct soundfield
        T(micCnt+1:micCnt+nbMic,:,J) = 1i.^(repmat(l.',nbMic,1)) ...
            .* ( a*Jkr(:,l+1,J).*Y - 1i*b*JkrDer(:,l+1,J).*Y.*Urad ...
            - 1i*b*Jkr(:,l+1,J)./repmat((k(J)+eps*(k(J)==0))*r.*cos(elv),1,nbComp).*repmat(m.',nbMic,1).*Y_.*Uazm ...
            + 1i*b*Jkr(:,l+1,J).*repmat(cos(elv)/(k(J)+eps*(k(J)==0))./r,1,nbComp).*convCoef.*Pder.*Truc.*Uelv ) ;
        
        if (freq(J)>0) && (R>0)
            
            % Admittance
            A = 1/(Z(freq(J))+eps*(Z(freq(J))==0)) ;
            
            % Values of the spherical Hankel functions
            % (Hl(kr) & Hl(kR))
            Hkr = SphericalBesselh( ...
                repmat(0:maxOrder,nbMic,1),2, ...
                repmat(k(J)*r,1,maxOrder+1) ) ;
            HkR = SphericalBesselh(0:maxOrder,2,k(J)*R) ;
            
            % Values of the spherical Hankel function derivatives
            % (Hl'(kr) & Hl'(kR))
            HkrDer = SphericalBesselhDerivative( ...
                repmat(0:maxOrder,nbMic,1),2, ...
                repmat(k(J)*r,1,maxOrder+1) ) ;
            HkRDer = SphericalBesselhDerivative(0:maxOrder,2,k(J)*R) ;
            
            % Values of jl(kR) and jl'(kR)
            JkR = SphericalBesselj(0:maxOrder,k(J)*R) ;
            JkRDer = SphericalBesseljDerivative(0:maxOrder,k(J)*R) ;
            
            % The matrix is corrected to take into account
            % the field diffracted by the sphere
            T(micCnt+1:micCnt+nbMic,:,J) = ...
                T(micCnt+1:micCnt+nbMic,:,J) ...
                - 1i.^(repmat(l.',nbMic,1)) ...
                .* ( a*Hkr(:,l+1).*Y ...
                - 1i*b*HkrDer(:,l+1).*Y.*Urad ...
                - 1i*b*Hkr(:,l+1)./repmat((k(J)+eps*(k(J)==0))*r.*cos(elv),1,nbComp).*repmat(m.',nbMic,1).*Y_.*Uazm ...
                + 1i*b*Hkr(:,l+1).*repmat(cos(elv)/(k(J)+eps*(k(J)==0))./r,1,nbComp).*convCoef.*Pder.*Truc.*Uelv ) ...
                .* repmat((JkRDer(l+1)+1i*A*JkR(l+1))./(HkRDer(l+1)+1i*A*HkR(l+1)),nbMic,1) ;
            
        end
        
        
    end
    
    % Apply the microphone frequency response if necessary
    if ~isempty(micFreqResp)
       T(micCnt+1:micCnt+nbMic,:,:) = T(micCnt+1:micCnt+nbMic,:,:) ...
           .* repmat(permute(micFreqResp(freq),[2 3 1]),[nbMic,nbComp]) ; 
    end
    
    % New micCnt value (current total number of microphones)
    micCnt = micCnt + nbMic ;
    
end

end