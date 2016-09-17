function hoa2SpkCfg = Hoa2SpkDecodingFilters(hoaFmt,spkFmt,hoa2SpkOpt)

if nargin < 3
    hoa2SpkOpt = GenerateHoa2SpkOpt ;
end

% Initialize the filter structure
if ~isempty(hoa2SpkOpt.sampFreq)
    nmbSmpFrq = length(hoa2SpkOpt.sampFreq) ;
    filters(nmbSmpFrq).sampFreq   = [] ;
    filters(nmbSmpFrq).nbInput    = [] ;
    filters(nmbSmpFrq).nbOutput   = [] ;
    filters(nmbSmpFrq).gainMatrix = [] ;
    filters(nmbSmpFrq).firMatrix  = [] ;
end

% Shortcuts
nmbHrm = hoaFmt.nbComp ;
hrmIdx = hoaFmt.index ;
nmbSpk = spkFmt.nbSpk ;
decTyp = hoa2SpkOpt.decodType ;
azmSpk = spkFmt.sphCoord(:,1) ;
elvSpk = spkFmt.sphCoord(:,2) ;
radSpk = spkFmt.sphCoord(:,3) ;
dstCor = hoa2SpkOpt.spkDistCor ;

% 2d/3d decoding
switch hoaFmt.conv
    case {'N2D','SN2D'}
        decDim = '2d' ;
    case {'N3D','SN3D'}
        decDim = '3d' ;
end


%%%%%%%%%%%%%%%%%%%%%%
% FILTER CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%

if dstCor
    
    % Create the filters for every sampling frequency value
    for I = 1 : nmbSmpFrq
        
        % Shortcuts
        smpFrq = hoa2SpkOpt.sampFreq(I) ;
        fltLng = hoa2SpkOpt.filterLength(I) ;
        trnFrq = hoa2SpkOpt.transFreq ;
        trnWdt = hoa2SpkOpt.transWidth ;
        
        % Number of frequency bins, frequency and wave number values
        fftNmb = 4*fltLng ;
        frq = [1/4,1:fftNmb/2]'*smpFrq/fftNmb ;
        wnb = 2*pi*frq/340 ;
        
        % Maximum HOA order value
        ord = max(hrmIdx(:,1)) ;
        
        % Spherical harmonic values for the speaker directions
        SphHrmSpk = SphericalHarmonicMatrix(hoaFmt,azmSpk,elvSpk) ;
        
        % "Modal coefficients" (radial components) for the
        % speakers, modelled as spherical source
        ModCoeSpk = Spk2HoaRadialWeightings(0:ord,radSpk*wnb') ;
        if nmbHrm == 1
            ModCoeSpk = permute(ModCoeSpk,[3 1 2]) ;
        end
        
        % Switch: mixed decoding / otherwise
        if strcmp(decTyp,'mixed')
            
            % Transition between basic and maxrE decoding
            if strcmp(trnFrq,'kr')
                
                % Transition frequency in Hz
                ord = max(hrmIdx(:,1)) ;
                trnFrq = 340/2/pi * (2*ord+1)/exp(1)/.1 + trnWdt/2 ;
                trnFrq = round(trnFrq) ;
                
                % Update the transition frequency the option structure
                hoa2SpkOpt.transFreq = trnFrq ;
                
            end
            trnIdx = round(trnFrq*fftNmb/smpFrq) + 1 ;
            trnLng = 2*round(trnWdt*fftNmb/smpFrq/2)+1 ;
            trnWdw = hann(2*trnLng-1) ;
            trnMax = [ zeros(trnIdx-(trnLng+1)/2,1) ; ...
                       trnWdw(1:trnLng) ; ...
                       ones(2*fltLng+1-trnIdx-(trnLng-1)/2,1) ] ;
            trnBas = 1 - trnMax ;
            
            % maxrE decoding matrix
            % If the speaker layout is quite unregular, then the matrix is
            % calculated so that it tries to match the e-VBAP panning
            ortErr = ...
                norm(eye(nmbHrm)-1/nmbSpk*(SphHrmSpk*SphHrmSpk.'),'fro') ...
                / sqrt(nmbHrm) ;
            if ortErr > .25
                if strcmp(decDim,'3d')
                    [azmVbp,elvVbp] = RegularPolyhedron(642) ;
                    Vbp = Vbap3dGainMatrix([azmSpk elvSpk],[azmVbp elvVbp], ...
                        'velocity','energy') ;
                else
                    azmVbp = linspace(0,358,180)'*pi/180 ;
                    elvVbp = zeros(180,1) ;
                    Vbp = Vbap2dGainMatrix(azmSpk,azmVbp, ...
                        'velocity','energy') ;
                end
                SphHrmVbp = SphericalHarmonicMatrix(hoaFmt,azmVbp,elvVbp) ;
                DecMax = Vbp * pinv(SphHrmVbp) ;
            else
                DecMax = Hoa2SpkDecodingMatrix(hoaFmt,azmSpk,elvSpk, ...
                    'decodType','maxrE','invMethod','tikhonov') ;
            end
            
            % Basic decoding matrix
            DecBas = Hoa2SpkDecodingMatrix(hoaFmt,azmSpk,elvSpk, ...
                'decodType','basic','invMethod','tikhonov') ;
            
            
            % Decoding filter frequency responses
            firMatrix = zeros(nmbSpk,nmbHrm,fftNmb) ;
            for J = 1 : fftNmb/2+1
                % Speaker distance correction matrix
                TrnMat = SphHrmSpk .* ModCoeSpk(hrmIdx(:,1)+1,:,J) ;
                TrnMatInv = TrnMat' * ((TrnMat*TrnMat' ...
                    +.01*norm(TrnMat,'fro')*eye(nmbHrm))^(-1)) ;
                SpkCor = TrnMatInv * SphHrmSpk ;
                % Mix the basic and maxrE decoding matrices and apply
                % the speaker distance correction
                firMatrix(:,:,J) = SpkCor ...
                    * (trnBas(J)*DecBas + trnMax(J)*DecMax)  ;
            end
            firMatrix = permute(firMatrix,[3 1 2]) ;
            firMatrix(fftNmb/2+2:end,:) = ...
                conj(firMatrix(fftNmb/2:-1:2,:)) ;
            
        else
            
            % basic / maxrE / inPhase weighting matrix
            if strcmp(decTyp,'maxrE') 
                WgtDec = Hoa2SpkMaxrEWeightings(ord,decDim) ;
            elseif strcmp(decTyp,'inPhase') 
                WgtDec = Hoa2SpkInPhaseWeightings(ord,decDim) ;
            else
                WgtDec = ones(nmbHrm,1) ;
            end
            WgtDec = sparse(1:nmbHrm,1:nmbHrm,WgtDec(hrmIdx(:,1)+1)) ;
            WgtDec = WgtDec * sqrt(nmbHrm) / norm(WgtDec,'fro') ;

            % Decoding filter frequency responses
            firMatrix = zeros(nmbSpk,nmbHrm,fftNmb) ;
            for J = 1 : fftNmb/2+1
                firMatrix(:,:,J) = pinv( SphHrmSpk ...
                    .* ModCoeSpk(hrmIdx(:,1)+1,:,J) ) * WgtDec  ;
            end
            firMatrix = permute(firMatrix,[3 1 2]) ;
            firMatrix(fftNmb/2+2:end,:) = ...
                conj(firMatrix(fftNmb/2:-1:2,:)) ;
                        
        end
        
        % Time-domain filters
        firMatrix = real(fftshift(ifft(firMatrix),1)) ;
        firMatrix = repmat(hanning(fltLng),[1 nmbSpk nmbHrm]) ...
            .* firMatrix(3*fftNmb/8+1:5*fftNmb/8,:,:) ;
        
        % hoa2MicFilters fields assignment
        filters(I) = struct('sampFreq',smpFrq, ...
            'nbInput',nmbHrm,'nbOutput',nmbSpk, ...
            'gainMatrix',[],'firMatrix',firMatrix) ;
        
    end
    
else
    
    % Switch: mixed decoding / otherwise
    if strcmp(decTyp,'mixed')
        
        % Creating the filters for every sampling ferquency value
        for I = 1 : nmbSmpFrq
            
            % Shortcuts
            smpFrq = hoa2SpkOpt.sampFreq(I) ;
            fltLng = hoa2SpkOpt.filterLength(I) ;
            trnFrq = hoa2SpkOpt.transFreq ;
            trnWdt = hoa2SpkOpt.transWidth ;
            
            % Number of frequency bins
            fftNmb = 4*fltLng ;
            
            % Basic decoding matrix
            DecBas = Hoa2SpkDecodingMatrix(hoaFmt, ...
                azmSpk,elvSpk,'decodType','basic') ;
            
            % maxrE decoding matrix
            % If the speaker layout is quite unregular, then the matrix is
            % calculated so that it tries to match the e-VBAP panning
            SphHrmSpk = SphericalHarmonicMatrix(hoaFmt,azmSpk,elvSpk) ;
            ortErr = ...
                norm(eye(nmbHrm)-1/nmbSpk*(SphHrmSpk*SphHrmSpk.'),'fro') ...
                / sqrt(nmbHrm) ;
            if ortErr > .25
                if strcmp(decDim,'3d')
                    [azmVbp,elvVbp] = RegularPolyhedron(642) ;
                    Vbp = Vbap3dGainMatrix([azmSpk elvSpk],[azmVbp elvVbp], ...
                        'velocity','energy') ;
                else
                    azmVbp = linspace(0,358,180)'*pi/180 ;
                    elvVbp = zeros(180,1) ;
                    Vbp = Vbap2dGainMatrix(azmSpk,azmVbp, ...
                        'velocity','energy') ;
                end
                SphHrmVbp = SphericalHarmonicMatrix(hoaFmt,azmVbp,elvVbp) ;
                DecMax = Vbp * pinv(SphHrmVbp) ;
            else
                DecMax = Hoa2SpkDecodingMatrix(hoaFmt,azmSpk,elvSpk, ...
                    'decodType','maxrE','invMethod','tikhonov') ;
            end
            
            % Transition between the decodings
            if strcmp(trnFrq,'kr')
                
                % Maximum HOA order
                ord = max(hrmIdx(:,1)) ;
                
                % Transition frequency in Hz
                trnFrq = 340/2/pi * (2*ord+1)/exp(1)/.1 + trnWdt/2 ;
                trnFrq = round(trnFrq) ;
                
                % Update the transition frequency the option structure
                hoa2SpkOpt.transFreq = trnFrq ;
                
            end
            trnIdx = round(trnFrq*fftNmb/smpFrq) + 1 ;
            trnLng = 2*round(trnWdt*fftNmb/smpFrq/2)+1 ;
            trnWdw = hann(2*trnLng-1) ;
            trnMax = [ zeros(trnIdx-(trnLng+1)/2,1) ; ...
                       trnWdw(1:trnLng) ; ...
                       ones(2*fltLng+1-trnIdx-(trnLng-1)/2,1) ] ;
            trnBas = 1 - trnMax ;
            
            % Playback filters frequency responses
            firMatrix = zeros(nmbSpk,nmbHrm,fftNmb) ;
            for J = 1 : fftNmb/2+1
                firMatrix(:,:,J) = ...
                    trnBas(J) * DecBas ...
                    + trnMax(J) * DecMax ;
            end
            firMatrix = permute(firMatrix,[3 1 2]) ;
            firMatrix(fftNmb/2+2:end,:) = ...
                conj(firMatrix(fftNmb/2:-1:2,:)) ;
            
            % Time-domain filters
            firMatrix = real(fftshift(ifft(firMatrix),1)) ;
            firMatrix = repmat(hanning(fltLng),[1 nmbSpk nmbHrm]) ...
                .* firMatrix(3*fftNmb/8+1:5*fftNmb/8,:,:) ;
            
            % hoa2MicFilters fields assignment
            filters(I) = struct( ...
                'sampFreq',hoa2SpkOpt.sampFreq(I), ...
                'nbInput',nmbHrm,'nbOutput',nmbSpk, ...
                'gainMatrix',[],'firMatrix',firMatrix) ;
            
        end
        
    else
        
        % Decoding matrix
        gainMatrix = Hoa2SpkDecodingMatrix(hoaFmt,azmSpk,elvSpk, ...
            'decodType',decTyp) ;
        
        % hoa2MicFilters fields assignment
        filters = struct('sampFreq',[], ...
            'nbInput',nmbHrm,'nbOutput',nmbSpk, ...
            'gainMatrix',gainMatrix,'firMatrix',[]) ;
        
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%
% FILTER EQUALIZATION %
%%%%%%%%%%%%%%%%%%%%%%%

% Equalise only when using FIR filters (as opposed to a decoding matrix)
if ~isempty(filters(1).firMatrix) 
    
    % Equalise the filters for every sampling frequency value
    for I = 1 : nmbSmpFrq
        
        % Shortcuts
        smpFrq = hoa2SpkOpt.sampFreq(I) ;
        fltLng = hoa2SpkOpt.filterLength(I) ;

        % Decoding filter FFTs and corresponding frequency bins
        DecFft = fft(fftshift(filters(I).firMatrix,1)) ;
        frq = [1/4,1:fltLng/2]'*smpFrq/fltLng ;
        
        % A set of plane-wave directions (2D/3D)
        if strcmp(decDim,'3d')
            [azmSou,elvSou] = RegularPolyhedron(642) ;
        else
            azmSou = (-180:2:178)'*pi/180 ;
            elvSou = zeros(size(azmSou)) ;
        end
        
        % order-L spherical harmonic expansion for these plane waves
        SphHrmExp = SphericalHarmonicMatrix(hoaFmt,azmSou,elvSou) ;
        
        % A set of points around the center of the loudspeaker array
        radPts = .12 ;
        if strcmp(decDim,'3d')
            stpPts = .02 ;
            xPts = (-radPts:stpPts:radPts)' ;
            yPts = xPts ;
            zPts = xPts ;
            [xPts,yPts,zPts] = meshgrid(xPts,yPts,zPts) ;
            xyzPts = [ xPts(:) yPts(:) zPts(:) ] ;
            xyzPts = xyzPts(sqrt(sum(xyzPts.^2,2))<=radPts,:) ;
        else
            stpPts = .01 ;
            xPts = (-radPts:stpPts:radPts)' ;
            yPts = xPts ;
            [xPts,yPts] = meshgrid(xPts,yPts) ;
            xyzPts = [ xPts(:) yPts(:) zeros(size(xPts(:))) ] ;
            xyzPts = xyzPts(sqrt(sum(xyzPts.^2,2))<=radPts,:) ;
        end
        
        % Equalise the filter so that the average energy of the
        % reconstructed sound field is the same for every frequency
        for J = 1 : fltLng/2+1
            
            % Transfers between the loudspeakers and the points
            % NB: the speaker gains and phases are supposed to be equalized
            SpkPts = SphericalWave(spkFmt.xyzCoord,xyzPts,frq(J),1) ;
        
            % Reconstructed sound field at the Jth frequency
            SndFld = SpkPts * permute(DecFft(J,:,:),[2 3 1]) * SphHrmExp ;
            
            % Filter equalization
            DecFft(J,:,:) = DecFft(J,:,:) / norm(SndFld,'fro') ; 
            
        end
        DecFft(fltLng/2+2:fltLng,:,:) = conj(DecFft(fltLng/2:-1:2,:,:)) ;

        % Back to the time domain
        filters(I).firMatrix = fftshift(real(ifft(DecFft)),1) ;
        
        % General normalisation: maximum amplification is 0dB
        filters(I).firMatrix = filters(I).firMatrix ...
            / max(max(max(abs(fft(filters(I).firMatrix))))) ;

    end

end


%%%%%%%%%%
% OUTPUT %
%%%%%%%%%%

% Fill the hoa2SpkCfg structure
hoa2SpkCfg.hoaFmt     = hoaFmt ;
hoa2SpkCfg.spkFmt     = spkFmt ;
hoa2SpkCfg.hoa2SpkOpt = hoa2SpkOpt ;
hoa2SpkCfg.filters    = filters ;