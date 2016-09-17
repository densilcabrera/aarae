function transfers = Spk2MicTransfers(spkFmt,micFmt,frq,radRef)

% Default radRef is empty
if nargin<4
    radRef = [] ;
end

% Format the frequency values.
frq = frq(:) ;

% Corresponding wave number values
wnb = 2*pi*frq/340 ;

% Number of sources and microphones
nmbSpk = spkFmt.nbSpk ;
nmbMic = micFmt.nbMic ;
nmbFrq = length(frq) ;

% Converting the 'reference radius' in complex speaker gains
if ~isempty(radRef)
    gaiSpk = radRef * exp(1i*wnb*(radRef*ones(1,nmbSpk))) ...
        .* exp(-1i*wnb*spkFmt.sphCoord(:,3)') ... 
        ./ repmat(spkFmt.sphCoord(:,3).',nmbFrq,1) ;
    gaiSpk(:,isinf(spkFmt.sphCoord(:,3))) = 1 ;
else
    gaiSpk = ones(nmbFrq,nmbSpk) ;
end


% if the microphone array is baffled or if the current sub-array consist in
% 'measured' mics, then the transfers have to be modelled using the
% spherical harmonic expansion. Otherwise, the PlaneWave and SphericalWave
% functions are used.
if (micFmt.sphereRadius>0)
    
    % Maximum hoa order required to model the microphone at fs/2
    hoaOrd = Hoa2MicRequiredOrder(micFmt,max(frq)) ;
    hoaFmt = GenerateHoaFmt('res2d',hoaOrd,'res3d',hoaOrd) ;
    nmbHrm = hoaFmt.nbComp ;
    
    % Spherical transfer matrix
    SphTrnMat = Hoa2MicTransferMatrix(hoaFmt,micFmt,frq) ;

    % Speaker spherical harmonic expansion coefficients
    OrdWgtSpk = Spk2HoaRadialWeightings(0:hoaOrd, ...
        wnb*spkFmt.sphCoord(:,3).') ;
    SphHarSpk = SphericalHarmonicMatrix(hoaFmt, ...
        spkFmt.sphCoord(:,1),spkFmt.sphCoord(:,2)) ;
    
    % Frequency responses
    % In the case where the sources are spherical waves at different
    % distances, these distances are taken into account in terms of gain
    % and delays
    transfers = zeros(nmbMic,nmbSpk,nmbFrq) ;
    for I = 1 : nmbFrq
        transfers(:,:,I) = SphTrnMat(:,:,I) ...
            * sparse(1:nmbHrm,1:nmbHrm,OrdWgtSpk(hoaFmt.index(:,1)+1,I)) ...
            * SphHarSpk * diag(gaiSpk(I,:)) ;
    end
    
else
    
    % microphone counter
    micCnt = 0 ;
    
    % Initialise the impulse response array
    transfers = zeros(nmbMic,nmbSpk,nmbFrq) ;
    
    for I  = 1 : micFmt.nbArrays
        
        % Current microphone index
        nmbMicCur = micFmt.arrays(I).nbMic ;
        micIdx = micCnt+1:micCnt + nmbMicCur ;
        
        % Current microphone positions (in the sub-array)
        xyzMic = micFmt.arrays(I).xyzCoord ;
        
        % Speaker position shortcut
        xyzSpk = spkFmt.xyzCoord ;
        
        if strcmp(micFmt.arrays(I).micType,'measured')
            
            % MicFmt corresponding to the current sub-array
            curMicFmt = GenerateMicFmt({'micType','measured', ...
                'measurements',micFmt.arrays(I).measurement, ...
                'xyzCoord',xyzMic}) ;
            
            % Maximum hoa order required to model the microphone at fs/2
            hoaOrd = Hoa2MicRequiredOrder(curMicFmt,smpFrq/2) ;
            hoaFmt = GenerateHoaFmt('res2d',hoaOrd,'res3d',hoaOrd) ;
            nmbHrm = hoaFmt.nbComp ;
            
            % Spherical transfer matrix
            SphTrnMat = Hoa2MicTransferMatrix(hoaFmt,curMicFmt,frq) ;
            
            % Speaker spherical harmonic expansion coefficients
            OrdWgtSpk = Spk2HoaRadialWeightings(0:hoaOrd,wnb*radSpk.') ;
            SphHarSpk = SphericalHarmonicMatrix(hoaFmt,azmSpk,elvSpk) ;
            
            % Frequency responses
            % In the case where the sources are spherical waves at
            % different distances, these distances are taken into account
            % in terms of gain and delays
            for J = 1 : nmbFrq
                transfers(micIdx,:,J) = SphTrnMat(:,:,J) ...
                    * sparse(1:nmbHrm,1:nmbHrm, ...
                    OrdWgtSpk(hoaFmt.index(:,1)+1,J)) ...
                    * SphHarSpk * diag(gaiSpk(J,:))  ;
            end
            
        else
            
            % If the current sub-array consists in 'omni' mics, the
            % transfers are computed using only PlaneWave or SphericalWave.
            % Otherwise, we need to take into account the microphone
            % directivity.
            
            % Transfer functions in the 'omni' case
            if ~any(isinf(spkFmt.sphCoord(:,3)))
                
                for J = 1 : nmbFrq
                    transfers(micIdx,:,J) = ...
                        SphericalWave(xyzSpk,xyzMic,frq(J),1) ...
                        * diag(gaiSpk(J,:)) ;
                end
                
            else
                [spkDir(:,1),spkDir(:,2),spkDir(:,3)] = ...
                    sph2cart(spkFmt.sphCoord(:,1),spkFmt.sphCoord(:,2),1) ;
                for J = 1 : nmbFrq
                    transfers(micIdx,:,J) = ...
                        PlaneWave(spkDir,xyzMic,frq(J))  ;
                end
                
            end
            
            % Directive microphone case
            if ~strcmp(micFmt.arrays(I).micType,'omni') ;
                
                % Spherical wave / Plane wave
                if ~any(isinf(spkFmt.sphCoord(:,3)))
                    
                    % Speaker position relative to every microphone
                    relXyzPos = repmat(permute(reshape(xyzSpk,nmbSpk,1,3), ...
                        [2 1 3]),nmbMicCur,1) ...
                        - repmat(reshape(xyzMic,nmbMicCur,1,3),1,nmbSpk) ;
                    
                    % Corresponding angles
                    [relSphPos(:,:,1),relSphPos(:,:,2)] = ...
                        cart2sph(relXyzPos(:,:,1),relXyzPos(:,:,2), ...
                        relXyzPos(:,:,3)) ;
                    
                else
                    
                    % Plane wave angles
                    relSphPos(:,:,1) =  ...
                        repmat(spkFmt.sphCoord(:,1).',nmbMicCur,1) ;
                    relSphPos(:,:,2) =  ...
                        repmat(spkFmt.sphCoord(:,2).',nmbMicCur,1) ;
                    
                end
                
                % Angles relative to the microphone directions
                micDir = micFmt.arrays(I).sphDirect ;
                relSphPos(:,:,1) = relSphPos(:,:,1) ...
                    - repmat(micDir(:,1),1,nmbSpk) ;
                relSphPos(:,:,2) = relSphPos(:,:,2) ...
                    - repmat(micDir(:,2),1,nmbSpk) ;
                
                % Microphone directivity function
                switch micFmt.arrays(I).micType
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
                DirFun = @(azm,elv) a + b * cos(azm).*cos(elv) ;
                
                % Applying the directivity
                transfers(micIdx,:,:) = transfers(micIdx,:,:) ...
                    .* repmat(DirFun(relSphPos(:,:,1),relSphPos(:,:,2)), ...
                    [1 1 nmbFrq]) ;
                
            end
            
        end
        
        % Updating the microphone counter
        micCnt = micCnt + nmbMicCur ;
        
    end
    
end

end
