function mic2HoaCfg = Mic2HoaEncodingFilters(hoaFmt,micFmt,mic2HoaOpt)

% Default encoding options
if nargin < 3
    mic2HoaOpt = GenerateMic2HoaOpt ;
end

% Number of sub-antennas
nmbArr = micFmt.nbArrays ;

% Test if the micFmt is compatible with post-Equalization if required.
% If not, the filtering method is switched to firMatrix.
if strcmp(mic2HoaOpt.filterType,'postEqFir')
    micTypes = cell(nmbArr,1) ;
    for I = 1 : nmbArr
        micTypes{I} = micFmt.arrays(I).micType ;
    end
    if ~any(~strcmp(micTypes,micTypes{1}))
        micType = micFmt.arrays(1).micType ;
        if strcmp(micType,'measured') && (nmbArr > 1)
                fprintf(['WARNING: ' ...
                    'Post-equalization cant be used for ' ...
                    'different arrays of measured mics... ' ...
                    'switching to FIR matrix method \n']) ;
                mic2HoaOpt.filterType = 'firMatrix' ;
                mic2HoaOpt.higherOrders  = false ;
                mic2HoaOpt.subArrayFilt  = false ;
        elseif strcmp(micType,'measured') && ...
            (size(micFmt.arrays.measurements.impulseResponses,3) > 1)
                fprintf(['WARNING: ' ...
                    'Post-equalization cant be used for ' ...
                    'an array of measured mics defined by a 3D ' ...
                    'impulse response array... ' ...
                    'switching to FIR matrix method \n']) ;
                mic2HoaOpt.filterType = 'firMatrix' ;
                mic2HoaOpt.higherOrders  = false ;
                mic2HoaOpt.subArrayFilt  = false ;
        else
            micDist = [] ;
            for I = 1 : nmbArr
                micDist = [micDist;micFmt.arrays(I).sphCoord(:,3)] ; %#ok<AGROW>
            end
            if ~any(abs(micDist-micDist(1))>1e-3)
                micPos = [] ;
                for I = 1 : nmbArr
                    micPos = [micPos;micFmt.arrays(I).sphCoord(:,1:2)] ; %#ok<AGROW>
                end
                micDir = [] ;
                for I = 1 : nmbArr
                    micDir = [micDir;micFmt.arrays(I).sphDirect(:,1:2)] ; %#ok<AGROW>
                end
                if any(any(micPos~=micDir))
                    fprintf(['WARNING: ' ...
                        'Microphones must be centrifugally ' ...
                        'oriented to use post-equalization...' ...
                        ' switching to FIR matrix method \n']) ;
                    mic2HoaOpt.filterType = 'firMatrix' ;
                    mic2HoaOpt.higherOrders  = false ;
                    mic2HoaOpt.subArrayFilt  = false ;
                end
            else
                fprintf(['WARNING: ' ...
                    'Microphones must all be the same distance ' ...
                    'from the origin to use post-equalization...' ...
                    ' switching to FIR matrix method \n']) ;
                mic2HoaOpt.filterType = 'firMatrix' ;
                mic2HoaOpt.higherOrders  = false ;
                mic2HoaOpt.subArrayFilt  = false ;
            end
        end
    else
        fprintf(['WARNING: ' ...
            'Microphones must all be the same type to use ' ...
            'post-equalization... switching to FIR matrix method \n']) ;
        mic2HoaOpt.filterType = 'firMatrix' ;
        mic2HoaOpt.higherOrders  = false ;
        mic2HoaOpt.subArrayFilt  = false ;
    end
end

% Test if sub-array filtering is appropriate
if mic2HoaOpt.subArrayFilt
    if nmbArr == 1
        fprintf(['WARNING: ' ...
            'There is no point using sub-array filtering with a ' ...
            'simple array... switching subArrayFilt to false \n']) ;
        mic2HoaOpt.subArrayFilt = false ;
    end
end
        
% filter structure initialization
nbSampFreq = length(mic2HoaOpt.sampFreq) ;
filters(nbSampFreq).sampFreq   = [] ;
filters(nbSampFreq).nbInput    = [] ;
filters(nbSampFreq).nbOutput   = [] ;
filters(nbSampFreq).gainMatrix = [] ;
filters(nbSampFreq).postEqFir  = [] ;
filters(nbSampFreq).firMatrix  = [] ;

% Shortcuts
nbHarm    = hoaFmt.nbComp ;
nbMic     = micFmt.nbMic ;
maxOrder  = max(hoaFmt.index(:,1)) ; 
lowPasFrq = mic2HoaOpt.lowPassFreq ;

% Generate a limiting function depending on the limiter method
switch mic2HoaOpt.limiterMethod
    case 'tikh'
        LimInv = @(M,limCoef) conj(M) ./ (abs(M).^2 + (1/2*limCoef)^2) ;        
    case 'expo'
        LimInv = @(M,limCoef) (1-exp(-(max(abs(M),eps)/(0.6382*limCoef)).^2)) ...
            .* exp(-1i*angle(M)) ./ max(abs(M),eps) ;
    case 'max'
        LimInv = @(M,limCoef) exp(-1i*angle(M)).*min(1./abs(M),1/limCoef) ;
end

% Intermediate "hoaFmt" used in the calculations
hoaFmtCal = GenerateHoaFmt('res2d',hoaFmt.res2d,'res3d',hoaFmt.res3d) ;

% Conversion matrix from hoaFmtCal to hoaFmt
CnvMat = Hoa2HoaConversionMatrix(hoaFmtCal,hoaFmt) ;

% Compute the filters depending on the filtering method
switch mic2HoaOpt.filterType
    
    case 'postEqFir'
        
        % Compute the matrix of spherical harmonic coefficients
        Y = SphericalHarmonicMatrix(hoaFmtCal,micPos(:,1),micPos(:,2)) ;
         
        % Compute the gain matrix
        gainMatrix = pinv(Y.') ;
        
        % Compute the post-equalization filters for the different sampling
        % frequencies
        for I = 1 : nbSampFreq
            
            % Useful shortcuts
            fftNb = mic2HoaOpt.filterLength(I) ;
            fs    = mic2HoaOpt.sampFreq(I) ;
            
            % Time-domain window
            timeWindow = hanning(fftNb) ;
            
            % Vector of the frequency values of interest
            freq  = (0:fs/fftNb:fs/2)' ;
            
            % Compute the radial weightings (modal coefficients)
            if strcmp(micType,'measured')
                W = Hoa2MicRadialWeightingsFromMeas(0:maxOrder,freq, ...
                    micFmt.arrays(1).sphCoord(1,3), ...
                    micFmt.arrays(1).measurements) ;
            else
                W = Hoa2MicRadialWeightings(0:maxOrder, ...
                    2*pi*freq/340*micDist(1),'micType',micType, ...
                    'kR',2*pi*freq/340*micFmt.sphereRadius, ...
                    'impedance',micFmt.sphereImpedance(freq), ...
                    'freqResp',micFmt.arrays(1).micFreqResp(freq)) ;
            end
            
            % Low-pass window to avoid oscillations at high frequencies
            nbPtsTrn = round((fs/2-lowPasFrq)/fs*fftNb) ;
            wnd = hanning(2*nbPtsTrn,'periodic') ;
            lowPassWindow = [ones(fftNb/2+1-nbPtsTrn,1);wnd(end/2+1:end)] ;
            
            % Compute the regularization coefficient
            limCoef = 1 / sqrt(micFmt.nbMic) ...
                / 10^(mic2HoaOpt.limiterLevel/20) ;
            
            % Encoding filters responses
            postEqFirFft = zeros(fftNb,size(W,1)) ;
            postEqFirFft(1:fftNb/2+1,:) = LimInv(W,limCoef).' ...
                .* repmat(lowPassWindow,1,size(W,1)) ;
            
            % High frequency equalization if required
            if mic2HoaOpt.highFreqEq
                % Complete matrix T (including high orders)
                ordPhy = Hoa2MicRequiredOrder(micFmt,max(freq)) ;
                hoaFmtPhy = GenerateHoaFmt('res2d',ordPhy,'res3d',ordPhy) ;
                T = Hoa2MicTransferMatrix(hoaFmtPhy,micFmt,freq) ;
                % Filter responses corresponding to the antenna behaviour
                % (acoustic diffraction/directivity + encoding filters)
                micBehavFft = ...
                    zeros(hoaFmt.nbComp,hoaFmtPhy.nbComp,length(freq)) ;
                for J = 1 : fftNb/2+1
                    micBehavFft(:,:,J) = ...
                        ( diag(postEqFirFft(J,hoaFmt.index(:,1)+1)) ...
                        * gainMatrix ) * T(:,:,J) ;
                end
                % List of the orders in hoaFmt
                for J = 0 : max(hoaFmt.index(:,1))
                    if any(hoaFmt.index(:,1)==J)
                        hoaOrder(J+1) = J ; %#ok<AGROW>
                    end
                end
                % Matrix of hoaFmt and hoaFmtPhy common harmonics
                commonHarm = repmat(permute(hoaFmt.index,[1 3 2]), ...
                    1,hoaFmtPhy.nbComp) ...
                    == repmat(permute(hoaFmtPhy.index,[3 1 2]), ...
                    hoaFmt.nbComp,1) ;
                commonHarm = commonHarm(:,:,1) & commonHarm(:,:,2) ;
                % Encoding error for every harmonic at every frequency
                encError = permute(sum(abs(micBehavFft ...
                    - repmat(commonHarm,[1 1 length(freq)])).^2,2),[1 3 2]) ;
                % Estimate the spatial aliasing frequency
                aliasFreqInd = ones(length(hoaOrder),1) ;
                aliasFreq    = freq(aliasFreqInd) ;
                for J = 1 : length(hoaOrder)
                    encErrJthOrd = ...
                        encError(hoaFmt.index(:,1)== hoaOrder(J),:) ;
                    if any(mean(encErrJthOrd,1)<=.1)
                        aliasFreqInd(J) = ...
                            find(mean(encErrJthOrd,1)<=.1,1,'last') ;
                        aliasFreq(J) = freq(aliasFreqInd(J)) ;
                    elseif J~=1
                        aliasFreq(J) = aliasFreq(J-1) ;
                        aliasFreqInd(J) = aliasFreqInd(J-1) ;
                    end
                end
                % Equalization gains
                equalGains = zeros(fftNb/2+1,length(hoaOrder)) ;
                for J = 1 : length(hoaOrder)
                    jthOrdHarms = hoaFmt.index(:,1)==hoaOrder(J) ;
                    equalGains(:,J) = sqrt(sum(sum(abs( ...
                        repmat(permute(lowPassWindow,[2 3 1]), ...
                        length(find(jthOrdHarms)),hoaFmtPhy.nbComp) ...
                        .* repmat(commonHarm(jthOrdHarms,:), ...
                        [1 1 fftNb/2+1])).^2,2),1) ...
                        ./ (sum(sum(abs(micBehavFft( ...
                        jthOrdHarms,:,:)).^2,2),1)+eps) ) ;
                end
                % Equalize the encoding filters
                lowPsFreqInd = round(lowPasFrq/fs*fftNb) ;
                lowPsTrnLng = min([2*round(1000/fs*fftNb), ...
                    2*floor(min(lowPsFreqInd)/2), ...
                    2*floor((fftNb/2+1-lowPsFreqInd))]) ;
                equalPostEqFir = zeros(fftNb/2+1,length(hoaOrder)) ;
                for J = 1 : length(hoaOrder)
                    aliasTrnLng = min([2*round(1000/fs*fftNb), ...
                        2*floor(min(aliasFreqInd)/2), ...
                        2*floor((fftNb/2+1-aliasFreqInd(J)))]) ;
                    equalCoef = [ zeros(aliasFreqInd(J)-aliasTrnLng/2,1) ; ...
                        sin(pi/2/(aliasTrnLng-1)*(0:aliasTrnLng-1)').^2 ; ...
                        ones(length(freq)-aliasFreqInd(J)-aliasTrnLng/2,1)] ;
                    equalCoef = equalCoef ...
                        .* [ ones(lowPsFreqInd-lowPsTrnLng/2,1) ; ...
                        cos(pi/2/(lowPsTrnLng-1)*(0:lowPsTrnLng-1)').^2 ; ...
                        zeros(fftNb/2+1-lowPsFreqInd-lowPsTrnLng/2,1)] ;
                    gain = 1+(equalGains(:,J)-1).*equalCoef ;
                    equalPostEqFir(:,J) = ...
                        [ postEqFirFft(1:fftNb/2,J).*gain(1:end-1) ; 0 ] ;
                end
                postEqFirFft(1:fftNb/2+1,:) = equalPostEqFir ;
            end
            
            % Generate the time-domain filters
            postEqFirFft(fftNb/2+2:end,:) = ...
                conj(postEqFirFft(fftNb/2:-1:2,:)) ;
            postEqFir = real(fftshift(ifft(postEqFirFft),1)) ;
            postEqFir = repmat(timeWindow,[1 size(postEqFir,2)]) ...
                .* postEqFir ;
            
            % Apply the format conversion to the gain matrix
            gainMatrix = CnvMat * gainMatrix ;
            
            % Assign the fields in the encoding filter structure
            filters(I) = struct( ...
                'sampFreq',mic2HoaOpt.sampFreq(I), ...
                'nbInput',micFmt.nbMic,'nbOutput',hoaFmt.nbComp,...
                'gainMatrix',gainMatrix,'postEqFir',postEqFir, ...
                'firMatrix',[]) ;
            
        end
        
    case 'firMatrix'
        
        % Compute the regularization coefficient
        limCoef = 1 / 10^(mic2HoaOpt.limiterLevel/20) ;
            
        % Compute the encoding filters for the different sampling
        % frequencies
        for I = 1 : length(mic2HoaOpt.sampFreq)
            
            % Useful shortcuts
            fftNb = mic2HoaOpt.filterLength(I) ;
            fs    = mic2HoaOpt.sampFreq(I) ;
            
            % Time-domain window
            timeWindow = hanning(fftNb) ;
            
            % Vector of the frequency values of interest
            freq  = [ 1 ; (1:fftNb/2)'*fs/fftNb ] ;
            
            % Low-pass window to avoid oscillations at high frequencies
            nbPtsTrn = round(2000/fs*fftNb) ;
            wnd = hanning(2*nbPtsTrn,'periodic') ;
            lowPassWindow = [ones(fftNb/2+1-nbPtsTrn,1);wnd(end/2+1:end)] ;
            
            % Maximum order of the harmonics taken into account
            if mic2HoaOpt.higherOrders
                ordMod = Hoa2MicRequiredOrder(micFmt,max(freq)) ;
                hoaFmtMod = GenerateHoaFmt('res2d',ordMod,'res3d',ordMod) ;
            else
                hoaFmtMod = hoaFmtCal ;
            end
            
            % Compute matrix "T"
            T = Hoa2MicTransferMatrix(hoaFmtMod,micFmt,freq) ;
            
            % Depending on the 'subArrayFilt' option,
            % compute "matrix T"'s regularized pseudo-inverse, OR
            % compute the filter responses for each sub-array and then
            % mix them to achieve the best performance at every freq
            firMatrixFft = zeros(fftNb,nbHarm,nbMic) ;
            if ~mic2HoaOpt.subArrayFilt
                invS = zeros(hoaFmtMod.nbComp,nbMic) ;
                for J = 1 : fftNb/2+1
                    [U,S,V] = svd(T(:,:,J)) ;
                    invS(1:min(size(S)),1:min(size(S))) = ...
                        diag(LimInv(diag(S(1:min(size(S)), ...
                        1:min(size(S)))),limCoef)) ;
                    invT = V * invS * U' ;
                    firMatrixFft(J,:,:) = permute( lowPassWindow(J) ...
                        * invT(1:nbHarm,:) , [3 1 2] ) ;
                end
            else
                % We need the complete matrix T (including high orders)
                ordPhy = Hoa2MicRequiredOrder(micFmt,max(freq)) ;
                hoaFmtPhy = GenerateHoaFmt('res2d',ordPhy,'res3d',ordPhy) ;
                Tphy = Hoa2MicTransferMatrix(hoaFmtPhy,micFmt,freq) ;
                % Assumed RMS measurement noise level, and power
                noiLvl = -30 ;
                noiPow = 10^(noiLvl/10) ;
                % Calculate the theoretical encoding error from each array
                micCnt = 0 ;
                EncErr = zeros(fftNb/2+1,nbHarm,nmbArr) ;
                for J = 1 : nmbArr
                    micIdx = micCnt+1 : micCnt+micFmt.arrays(J).nbMic ;
                    for K = 1 : fftNb/2+1
                        % Theoretical encoding filter response
                        EncFltThe = pinv(T(micIdx,:,K)) ;
                        % Product of the Enc. filter matrix with Tphy 
                        MatPrd = EncFltThe * Tphy(micIdx,:,K) ;
                        for L = 1 : nbHarm
                            % Measurement noise error
                            NoiErr = sum(abs(EncFltThe(L,:)).^2) * noiPow ; 
                            % Aliasing error
                            AlsErr = sum(abs(MatPrd(L,nbHarm+1:end)).^2) ;
                            % Total encoding error 
                            EncErr(K,L,J) = NoiErr + AlsErr ;
                        end
                    end
                    micCnt = micCnt + micFmt.arrays(J).nbMic ;
                end
                % Calculate the optimal weights for the sub-arrays
                % (for every frequency and every spherical harmonic)
                ArrWgt = 1 ./ EncErr ;
                ArrWgt = ArrWgt ./ repmat(sum(ArrWgt,3),[1 1 nmbArr]) ;
                % Individual encoding filter responses for each sub-array
                % (Regularised/limited pseudo-inverse)
                micCnt = 0 ;
                for J = 1 : nmbArr
                    invS = zeros(hoaFmtMod.nbComp,micFmt.arrays(J).nbMic) ;
                    micIdx = micCnt+1 : micCnt+micFmt.arrays(J).nbMic ;
                    for K = 1 : fftNb/2+1
                        [U,S,V] = svd(T(micIdx,:,K)) ;
                        invS(1:min(size(S)),1:min(size(S))) = ...
                            diag(LimInv(diag(S(1:min(size(S)), ...
                            1:min(size(S)))),limCoef)) ;
                        invT = V * invS * U' ;
                        firMatrixFft(K,:,micIdx) = ...
                            permute( invT(1:nbHarm,:) , [3 1 2] ) ;
                    end
                    micCnt = micCnt + micFmt.arrays(J).nbMic ;
                end
                % Total array filter responses
                micCnt = 0 ;
                for J = 1 : nmbArr
                    micIdx = micCnt+1 : micCnt+micFmt.arrays(J).nbMic ;
                    firMatrixFft(1:fftNb/2+1,:,micIdx) = ...
                        repmat(ArrWgt(:,:,J), ...
                        [1 1 micFmt.arrays(J).nbMic]) ...
                        .* firMatrixFft(1:fftNb/2+1,:,micIdx) ;
                    micCnt = micCnt + micFmt.arrays(J).nbMic ;
                end
                
            end
            
            % High frequency equalization if required
            if mic2HoaOpt.highFreqEq
                if ~mic2HoaOpt.subArrayFilt
                    % then we need the complete T matrix, unless it has
                    % already been computed
                    if ~mic2HoaOpt.higherOrders
                        ordPhy = Hoa2MicRequiredOrder(micFmt,max(freq)) ;
                        hoaFmtPhy = ...
                            GenerateHoaFmt('res2d',ordPhy,'res3d',ordPhy) ;
                        T = Hoa2MicTransferMatrix(hoaFmtPhy,micFmt,freq) ;
                    else
                        hoaFmtPhy = hoaFmtMod ; 
                    end
                else
                    % In the sub-array case,
                    % most things have already been computed
                    T = Tphy ;
                end
                % We also need the matrix of the common harmonics
                % between hoaFmt and hoaFmtPhy
                commonHarm = repmat(permute(hoaFmt.index,[1 3 2]), ...
                    1,hoaFmtPhy.nbComp) ...
                    == repmat(permute(hoaFmtPhy.index,[3 1 2]), ...
                    hoaFmt.nbComp,1) ;
                commonHarm = commonHarm(:,:,1) & commonHarm(:,:,2) ;
                % List of the orders in hoaFmt
                for J = 0 : max(hoaFmt.index(:,1))
                    if any(hoaFmt.index(:,1)==J)
                        hoaOrder(J+1) = J ; %#ok<AGROW>
                    end
                end
                % Filter responses corresponding to the antenna behaviour
                % (acoustic diffraction/directivity + encoding filters)
                micBehavFft = ...
                    zeros(hoaFmt.nbComp,hoaFmtPhy.nbComp,length(freq)) ;
                for J = 1 : fftNb/2+1
                    micBehavFft(:,:,J) = ...
                        permute(firMatrixFft(J,:,:),[2 3 1]) * T(:,:,J) ;
                end
                % Encoding error for every harmonic at every frequency
                encError = permute(sum(abs(micBehavFft ...
                    - repmat(commonHarm,[1 1 length(freq)])).^2,2),[1 3 2]) ;
                % For each order, find the spatial aliasing frequency
                aliasFreqInd = ones(length(hoaOrder),1) ;
                aliasFreq    = freq(aliasFreqInd) ;
                for J = 1 : length(hoaOrder)
                    encErrJthOrd = ...
                        encError(hoaFmt.index(:,1)== hoaOrder(J),:) ;
                    if any(mean(encErrJthOrd,1)<=.1)
                        aliasFreqInd(J) = ...
                            find(mean(encErrJthOrd,1)<=.1,1,'last') ;
                        aliasFreq(J) = freq(aliasFreqInd(J)) ;
                    elseif J~=1
                        aliasFreq(J) = aliasFreq(J-1) ;
                        aliasFreqInd(J) = aliasFreqInd(J-1) ;
                    end
                end
                % Equalization gains
                equalGains = zeros(fftNb/2+1,length(hoaOrder)) ;
                for J = 1 : length(hoaOrder)
                    jthOrdHarms = hoaFmt.index(:,1)==hoaOrder(J) ;
                    equalGains(:,J) = sqrt(sum(sum(abs( ...
                        repmat(permute(lowPassWindow,[2 3 1]), ...
                        length(find(jthOrdHarms)),hoaFmtPhy.nbComp) ...
                        .* repmat(commonHarm(jthOrdHarms,:), ...
                        [1 1 length(freq)])).^2,2),1) ...
                        ./ (sum(sum(abs(micBehavFft( ...
                        jthOrdHarms,:,:)).^2,2),1)+eps) ) ;
                end
                % Equalize the encoding filters
                lowPsFreqInd = round(lowPasFrq/fs*fftNb) ;
                lowPsTrnLng = min([2*round(1000/fs*fftNb), ...
                    2*floor(min(lowPsFreqInd)/2), ...
                    2*floor((fftNb/2+1-lowPsFreqInd))]) ;
                equalFirMatrix = zeros(fftNb/2+1,nbHarm,nbMic) ;
                for J = 1 : length(hoaOrder)
                    aliasTrnLng = min([2*round(1000/fs*fftNb), ...
                        2*floor(min(aliasFreqInd)/2), ...
                        2*floor((fftNb/2+1-aliasFreqInd(J)))]) ;
                    jthOrdHarms = hoaFmt.index(:,1)==hoaOrder(J) ;
                    equalCoef = [ zeros(aliasFreqInd(J)-aliasTrnLng/2,1) ; ...
                        sin(pi/2/(aliasTrnLng-1)*(0:aliasTrnLng-1)').^2 ; ...
                        ones(fftNb/2+1-aliasFreqInd(J)-aliasTrnLng/2,1) ] ;
                    equalCoef = equalCoef ...
                        .* [ ones(lowPsFreqInd-lowPsTrnLng/2,1) ; ...
                        cos(pi/2/(lowPsTrnLng-1)*(0:lowPsTrnLng-1)').^2 ; ...
                        zeros(fftNb/2+1-lowPsFreqInd-lowPsTrnLng/2,1)] ;
                    gain = 1 + (equalGains(:,J)-1) .* equalCoef ;
                    for K = 1 : nbMic
                        equalFirMatrix(:,jthOrdHarms,K) = ...
                            [ firMatrixFft(1:fftNb/2,jthOrdHarms,K) ...
                            .* repmat(gain(1:fftNb/2), ...
                            1,length(find(jthOrdHarms))) ; ...
                            zeros(1,length(find(jthOrdHarms))) ] ;
                    end
                end
                firMatrixFft(1:fftNb/2+1,:,:) = equalFirMatrix ;

            end
            
            % Convert from hoaFmtCal to hoaFmt
            for J = 1 : fftNb/2+1
                firMatrixFft(J,:,:) = ...
                    permute(CnvMat*squeeze(firMatrixFft(J,:,:)),[3 1 2]) ;
            end
            
            % Generate the time-domain filters
            firMatrixFft(fftNb/2+2:end,:,:) = ...
                conj(firMatrixFft(fftNb/2:-1:2,:,:)) ;
            firMatrix = fftshift(real(ifft(firMatrixFft)),1) ;
            firMatrix = firMatrix .* repmat(timeWindow, ...
                [1 size(firMatrixFft,2) size(firMatrixFft,3)]) ;
            
            % Assign the fields in the encoding filter structure
            filters(I) = struct( ...
                'sampFreq',mic2HoaOpt.sampFreq(I), ...
                'nbInput',micFmt.nbMic,'nbOutput',hoaFmt.nbComp,...
                'gainMatrix',[],'postEqFir',[],'firMatrix',firMatrix) ;
            
        end
        
end


% Output structure
mic2HoaCfg.hoaFmt     = hoaFmt ;
mic2HoaCfg.micFmt     = micFmt ;
mic2HoaCfg.mic2HoaOpt = mic2HoaOpt ;
mic2HoaCfg.filters    = filters ;
