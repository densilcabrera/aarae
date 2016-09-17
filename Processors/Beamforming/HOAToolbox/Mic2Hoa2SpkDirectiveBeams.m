function [achBeams,refBeams] = ...
    Mic2Hoa2SpkDirectiveBeams(mic2HoaCfg,Dec,frq,varargin)

% Default option values
noiseLvl = [] ;
plotOption  = true ;
dimMode     = '3d' ;
polarScale  = 'log' ;
angularRes  = 2 ;
souSphCoordPassed = false ;

% Number of optional arguments is checked
if round(length(varargin)/2)~=length(varargin)/2
    error('illegal number of arguments') ;
end

% New option values are checked and assigned
for I = 1 : 2 : length(varargin)-1
    switch varargin{I}
        case 'noiseLvl'
            if isreal(varargin{I+1}) && isscalar(varargin{I+1})
                noiseLvl = varargin{I+1} ;
            else
                error('noiseLvl must be a real value') ;
            end
        case 'plotOption'
            if islogical(varargin{I+1})
                plotOption = varargin{I+1} ;
            else
                error('plotOption must be a logical value') ;
            end
        case 'polarScale'
            switch varargin{I+1}
                case {'linear','log'}
                    polarScale = varargin{I+1} ;
                otherwise
                    error(['Polar plot scale available options ' ...
                        'are ''linear'' and ''log''']) ;
            end
        case 'dimMode'
            if strcmp(varargin{I+1},'2d') || strcmp(varargin{I+1},'3d')
                dimMode = varargin{I+1} ;
            else
                error('Unknown plotting mode') ;
            end
        case 'angularRes'
            if isreal(varargin{I+1}) && isscalar(varargin{I+1})
                angularRes = varargin{I+1} ;
            else
                error('Angular resolution must be a real value') ;
            end
        case 'souSphCoord'
            if isreal(varargin{I+1}) && (size(varargin{I+1},2)==2)
                souSphCoord = varargin{I+1} ;
            else
                error(['Source spherical coordinates must be passed ' ...
                    'in a Nx2 array of real values']) ;
            end
        otherwise
            error('Unknown option') ;
    end
end

% Shortcuts to the different parts of the mic2HoaCfg
hoaFmtEnc  = mic2HoaCfg.hoaFmt ;
micFmt     = mic2HoaCfg.micFmt ;
mic2HoaOpt = mic2HoaCfg.mic2HoaOpt ;
encFlt     = mic2HoaCfg.filters ;

% The encoding error is only computed for the highest sampling frequency
% value available
[smpFrqEnc,idxSmpFrq] = max(mic2HoaOpt.sampFreq) ;
fftNmb = mic2HoaOpt.filterLength(idxSmpFrq) ;
encFlt = encFlt(idxSmpFrq) ;

% Just a little test
if max(frq) > smpFrqEnc/2
    error('f cannot be higher than half the sampling frequency') ;
end

% HoaFmt required for an accurate simulation
ordPhy = Hoa2MicRequiredOrder(micFmt,max(frq)) ;
ordPhy = max(max(hoaFmtEnc.index(:,1)),ordPhy) ;
hoaFmtPhy = GenerateHoaFmt('res2d',ordPhy,'res3d',ordPhy, ...
    'type',hoaFmtEnc.type,'conv',hoaFmtEnc.conv) ;

% hoaFmtEnc -> hoaFmtPhys conversion matrix
PhyEncCnvMat = Hoa2HoaConversionMatrix(hoaFmtPhy,hoaFmtEnc) ;

% Encoding filter frequency responses
switch mic2HoaOpt.filterType
    case 'postEqFir'
        % Equalization filter Fourier Transforms are computed
        postEqFir = fft(fftshift(encFlt.postEqFir,1),[],1) ;
        postEqFir = postEqFir(:,hoaFmtEnc.index(:,1)+1) ;
        % Global encoding filter responses (including gain matrix)
        encFltRsp = repmat(postEqFir,[1 1 micFmt.nbMic]) ...
            .* repmat(permute(encFlt.gainMatrix,[3 1 2]),[fftNmb 1 1]) ;
    case 'firMatrix'
        % Encoding filter Fourier Transforms are computed
        encFltRsp = fft(fftshift(encFlt.firMatrix,1),[],1) ;
end

% Encoding filter responses at every 'frq' value
fftFrq = (0:fftNmb-1).'*smpFrqEnc/fftNmb ;
encFltRsp = interp1(fftFrq,encFltRsp,frq) ;
encFltRsp = permute(encFltRsp,[2 3 1]) ;

% Check the decoding argument, and compute the value of the "decoding"
% matrix for every frequency value
if isstruct(Dec)
    % Then it must be a hoa2Spk structure...
    % Decoding hoaFmt
    hoaFmtDec = Dec.hoaFmt ;
    % Decoding matrix/filters
    switch Dec.hoa2SpkOpt.decodType
        case 'mixed'
            % Then we have a matrix of filters
            % As for the encoding filters, we keep only the decoding
            % filters that have highest sampling frequency
            [smpFrqDec,idxSmpFrq] = max(Dec.hoa2SpkOpt.sampFreq) ;
            decFlt = Dec.filters(idxSmpFrq) ;
            % The response of the decoding filters 
            % is computed for every frequency
            fftNmb = Dec.hoa2SpkOpt.filterLength(idxSmpFrq) ;
            fftFrq = (0:fftNmb-1).'*smpFrqDec/fftNmb ;
            Dec = fft(fftshift(decFlt.firMatrix,1),[],1) ;
            Dec = interp1(fftFrq,Dec,frq) ;
            Dec = permute(Dec,[2 3 1]) ;
        otherwise
            % We have a decoding matrix
            Dec = Dec.filters.gainMatrix ;
            Dec = repmat(Dec,[1 1 length(frq)]) ;
    end
elseif isvector(Dec)||(length(size(Dec))<=2 && all(size(Dec)>1))
    % Then it should be a decoding "matrix" (each row represents a
    % speaker/beam direction and each column represents a sph. harmonic)
    hoaFmtDec = hoaFmtEnc ;
    if size(Dec,2)~=hoaFmtDec.nbComp 
        error(['The number of columns in the decoding matrix ' ...
            'doesn''t match the number of spherical harmonics'])
    end
    Dec = repmat(Dec,[1 1 length(frq)]) ;
elseif isempty(Dec)
    % Then a default set of weights is used, which corresponds to a "basic"
    % beam in the (0,0) direction
    hoaFmtDec = hoaFmtEnc ;
    Dec = SphericalHarmonicMatrix(hoaFmtDec,0,0).' ;
    Dec = repmat(Dec,[1 1 length(frq)]) ;
else
    % Sorry, but the code doesn't know what to do with that
    error('Unidentified decoding argument') ;    
end

% Conversion matrix between hoaFmtEnc and hoaFmtDec
EncDecCnvMat = Hoa2HoaConversionMatrix(hoaFmtEnc,hoaFmtDec) ;

% Source positions
if ~souSphCoordPassed
    if strcmp(dimMode,'2d')
        souSphCoord = (-180:angularRes:180).'*pi/180 ;
        souSphCoord = [souSphCoord zeros(size(souSphCoord))] ;
    else
        souAzm = (-180:angularRes:180)'*pi/180 ;
        souElv = (-090:angularRes:090)'*pi/180 ;
        [souAzm,souElv] = meshgrid(souAzm,souElv) ;
        souSphCoord = [ souAzm(:) souElv(:) ] ;
    end
elseif plotOption
    fprintf(['WARNING: source spherical coordinates have ' ...
        'been passed... plotOption is set to false\n']) ;
    plotOption = false ;
end
NmbSou = size(souSphCoord,1) ;

% Spherical harmonic components for the source and plot azimuth values
YSou = SphericalHarmonicMatrix(hoaFmtPhy, ...
    souSphCoord(:,1),souSphCoord(:,2)) ;

% HOA to mic transfer matrix
T = Hoa2MicTransferMatrix(hoaFmtPhy,micFmt,frq) ;

% Directivity Beams are computed, including measurement noise if required
nmbSpk = size(Dec,1) ;
achBeams = zeros(NmbSou,nmbSpk,length(frq)) ;
if ~isempty(noiseLvl)
    noi = rand(micFmt.nbMic,1) .* exp(1i*rand(micFmt.nbMic,1)*2*pi) ;
    noi = 10^(noiseLvl/20) * noi / sqrt(mean(abs(noi).^2)) ;
    N = diag(1+noi) ;    
    for I = 1 : length(frq)
        achBeams(:,:,I) = ( Dec(:,:,I) * EncDecCnvMat ...
            * encFltRsp(:,:,I) * N * T(:,:,I) * YSou ).'  ;
    end
else
    for I = 1 : length(frq)
        achBeams(:,:,I) = ( Dec(:,:,I) * EncDecCnvMat ...
            * encFltRsp(:,:,I) * T(:,:,I) * YSou ).'  ;
    end
end
refBeams = ( Dec(:,:,I) * EncDecCnvMat * PhyEncCnvMat * YSou ).' ;


% Plot the directivity beams
if plotOption
    
    for I = 1 : length(frq)
        
        % Maximum amplitude of the beams
        refBeamMax = max(max(abs(refBeams))) ;
        achBeamMax = max(max(abs(achBeams(:,:,I)))) ;
        
        for J = 1 : nmbSpk
            
            % Depending if we're in 2d or 3d, ...
            if strcmp(dimMode,'2d')
                
                circlesPlotAzim = (-pi : 2*pi/360 : pi)' ;
                
                % Polar plots are either using a log or a linear scale
                switch polarScale
                    case 'log'
                        ref = 20*log10(abs(refBeams(:,J)/refBeamMax))+30 ;
                        ach = 20*log10(abs(achBeams(:,J,I)/achBeamMax))+30 ;
                        ref(ref<0) = 0 ;
                        ach(ach<0) = 0 ;
                        maxLevel = 10*ceil(max(max(ref),max(ach))/10+1/10) ;
                        figure ;
                        hold on ;
                        for K = maxLevel : -10 : 10
                            patch(K*cos(circlesPlotAzim), ...
                                K*sin(circlesPlotAzim),[1 1 1], ...
                                'edgeColor',[.6 .6 .6],'linestyle',':');
                            text(-K*sqrt(2)/2,K*sqrt(2)/2, ...
                                [num2str(K-30) 'dB'], ...
                                'HorizontalAlignment','center', ...
                                'BackgroundColor',[1 1 1],...
                                'color',[.4 .4 .4]) ;
                        end
                    case 'linear'
                        ref = abs(refBeams(:,J)/refBeamMax) ;
                        ach = abs(achBeams(:,J,I)/refBeamMax) ;
                        maxLevel = ceil(max(max(ref),max(ach))*3+1/3)/3 ;
                        figure ;
                        hold on ;
                        for K = maxLevel : -1/3 : 1/3
                            patch(K*cos(circlesPlotAzim), ...
                                K*sin(circlesPlotAzim),...
                                [1 1 1],'edgeColor',[.6 .6 .6], ...
                                'linestyle',':');
                            switch round(round(K*3)/3) == round(K*3)/3
                                case 1
                                    text(-K*sqrt(2)/2,K*sqrt(2)/2, ...
                                        num2str(round(round(K*3)/3)), ...
                                        'HorizontalAlignment','center', ...
                                        'BackgroundColor',[1 1 1],...
                                        'color',[.4 .4 .4]) ;
                                case 0
                                    text(-K*sqrt(2)/2,K*sqrt(2)/2, ...
                                        [num2str(round(K*3)) '/3'], ...
                                        'HorizontalAlignment','center', ...
                                        'BackgroundColor',[1 1 1],...
                                        'color',[.4 .4 .4]) ;
                            end
                        end
                end
                line([0 0],[-maxLevel maxLevel],'color',[.9 .9 .9]) ;
                line([-maxLevel maxLevel],[0 0],'color',[.9 .9 .9]) ;
                plot(ref.*cos(souSphCoord(:,1)), ...
                    ref.*sin(souSphCoord(:,1)),'--b','linewidth',2) ;
                plot(ach.*cos(souSphCoord(:,1)), ...
                    ach.*sin(souSphCoord(:,1)),'r','linewidth',2) ;
                axis equal off ;
                hold off ;
                
            else
                
                % 3d plot
                ref = reshape(abs(refBeams(:,J)),size(souAzm)) ...
                    / refBeamMax ;
                ach = reshape(abs(achBeams(:,J,I)),size(souAzm)) ...
                    / refBeamMax ;
                [souX,souY,souZ] = sph2cart(souAzm,souElv,1) ;
                figure
                surf(ref.*souX,ref.*souY,ref.*souZ, ...
                    'EdgeColor',1/2*[0 0 1],'Facecolor',[0 0 1]) ;
                alpha(.5)
                hold on
                surf(ach.*souX,ach.*souY,ach.*souZ, ...
                    'EdgeColor',1/2*[1 .5 0],'FaceColor',[1 .5 0]) ;
                alpha(.5)
                axis equal ;
                light ;
                lighting phong ;
                campos([4 -4 1]*max(refBeamMax,achBeamMax))
                camproj perspective;
                
            end
            
        end
        
    end
    
end