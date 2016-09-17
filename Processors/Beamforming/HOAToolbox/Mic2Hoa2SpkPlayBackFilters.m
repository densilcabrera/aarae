function plBackFlt = Mic2Hoa2SpkPlayBackFilters(varargin)

%---------------------------%
% TESTING THE ARGUMENTS ETC %
%---------------------------%

% Test on the number of optional arguments
if round(size(varargin,2)/2)~=size(varargin,2)/2
    error('illegal number of arguments') ;
end

% Initialise arguments
micFmt    = [] ;
hoaFmtEnc = [] ;
encOpt    = [] ;
encFlt    = [] ;
hoaFmtDec = [] ;
spkFmt    = [] ;
decOpt    = [] ;
decFlt    = [] ;

% Initialise the flags
micFmtFlag     = 0 ;
hoaFmtEncFlag  = 0 ;
hoaFmtDecFlag  = 0 ;
spkFmtFlag     = 0 ;
encOptFlag     = 0 ;
decOptFlag     = 0 ;
mic2HoaCfgFlag = 0 ;
hoa2SpkCfgFlag = 0 ;

% Test and assign the optional arguments
for I = 1 : 2 : length(varargin)-1
    switch varargin{I}
        case 'micFmt'
            if isstruct(varargin{I+1})
                micFmt = varargin{I+1} ;
                micFmtFlag = 1 ;
            else
                error(['The microphone configuration must be passed in' ...
                    ' a structure. Please use the GenerateMicFmt' ...
                    'function']) ;
            end
        case 'hoaFmtEnc'
            if isstruct(varargin{I+1})
                hoaFmtEnc = varargin{I+1} ;
                hoaFmtEncFlag = 1 ;
            else
                error(['A set of spherical harmonics must be passed ' ...
                    'in a structure. Please use the GenerateHoaFmt' ...
                    'function']) ;
            end
        case 'hoaFmtDec'
            if isstruct(varargin{I+1})
                hoaFmtDec = varargin{I+1} ;
                hoaFmtDecFlag = 1 ;
            else
                error(['A set of spherical harmonics must be passed ' ...
                    'in a structure. Please use the GenerateHoaFmt' ...
                    'function']) ;
            end
        case 'spkFmt'
            if isstruct(varargin{I+1})
                spkFmt = varargin{I+1} ;
                spkFmtFlag = 1 ;
            else
                error(['The loudspeaker configuration must be passed ' ...
                    'in a structure. Please use the GenerateMic2HoaOpt' ...
                    'function']) ;
            end
        case 'encOpt'
            if isstruct(varargin{I+1})
                encOpt = varargin{I+1} ;
                encOptFlag = 1 ;
            else
                error(['Encoding options must be passed in a ' ...
                    'structure. Please use the GenerateMic2HoaOpt' ...
                    'function']) ;
            end
        case 'decOpt'
            if isstruct(varargin{I+1})
                decOpt = varargin{I+1} ;
                decOptFlag = 1 ;
            else
                error(['Decoding options must be passed in a ' ...
                    'structure. Please use the GenerateHoa2SpkOpt' ...
                    'function']) ;
            end
        case 'mic2HoaCfg'
            if isstruct(varargin{I+1})
                mic2HoaCfg = varargin{I+1} ;
                mic2HoaCfgFlag = 1 ;
            else
                error(['Encoding filters must be passed in a ' ...
                    'structure. Please use the Mic2HoaEncodingFilters' ...
                    'function']) ;
            end
        case 'hoa2SpkCfg'
            if isstruct(varargin{I+1})
                hoa2SpkCfg = varargin{I+1} ;
                hoa2SpkCfgFlag = 1 ;
            else
                error(['Decoding filters must be passed in a ' ...
                    'structure. Please use the Hoa2SpkDecodingFilters' ...
                    'function']) ;
            end
        otherwise
            error([varargin{I} ' is not a valid option']) ;
            
    end
end

% If a complete Mic -> HOA configuration has been passed, set the
% microphone array and encoding filters accordingly
if mic2HoaCfgFlag 
    if micFmtFlag || hoaFmtEncFlag || encOptFlag
        error(['A ''mic2HoaCfg'' cannot be passed at the same time as ' ...
            'a ''micFmt'', a ''hoaFmtEnc'' or a ''mic2HoaOpt'' structure.']) ;
    else
        micFmt    = mic2HoaCfg.micFmt ;
        hoaFmtEnc = mic2HoaCfg.hoaFmt ;
        encOpt    = mic2HoaCfg.mic2HoaOpt ;
        encFlt    = mic2HoaCfg.filters ;
    end
end

% If a complete HOA -> Spk configuration has been passed, set the
% loudspeaker array and decoding filters accordingly
if hoa2SpkCfgFlag 
    if hoaFmtDecFlag || spkFmtFlag || decOptFlag
        error(['A ''hoa2SpkCfg'' cannot be passed at the same time as ' ...
            'a ''hoaFmtDec'', a ''spkFmt'' or a ''hoaSpkOpt'' structure.']) ;
    else
        hoaFmtDec = hoa2SpkCfg.hoaFmt ;
        spkFmt    = hoa2SpkCfg.spkFmt ;
        decOpt    = hoa2SpkCfg.hoa2SpkOpt ;
        decFlt    = hoa2SpkCfg.filters ;
    end
end

% Microphone configuration if required
if isempty(micFmt)
    micFmt = GenerateMicFmt ;
end

% Encoding hoaFmt if required
if isempty(hoaFmtEnc)
    hoaFmtEnc = GenerateHoaFmt ;
end

% Encoding options if required
if isempty(encOpt)
    encOpt = GenerateMic2HoaOpt ;
end


% Encoding filters if required
if isempty(encFlt)
    encFlt = Mic2HoaEncodingFilters(hoaFmtEnc,micFmt,encOpt) ;
    encFlt = encFlt.filters ;
end

% Decoding hoaFmt if required
if isempty(hoaFmtDec)
    hoaFmtDec = hoaFmtEnc ;
end

% Loudspeaker configuration if required
if isempty(spkFmt)
    spkFmt = GenerateSpkFmt('nbSpk',2*max(hoaFmtDec.index(:,1))+2) ;
end

% Decoding options if required
if isempty(decOpt)
    if ~mic2HoaCfgFlag 
        % Encoding data structure
        mic2HoaCfg.micFmt     = micFmt ;
        mic2HoaCfg.hoaFmt     = hoaFmtEnc ;
        mic2HoaCfg.mic2HoaOpt = encOpt ;
        mic2HoaCfg.filters    = encFlt ;
    end
    % Evaluate a frequency for the basic to maxrE transition,
    % based on the expected encoding quality
    [EncSnr,snrFrq] = Mic2HoaSignalToNoiseRatio(mic2HoaCfg, ...
        'mode','orders','plotOption',false) ;
    EncOrd = repmat(0:size(EncSnr,2)-1,length(snrFrq),1) ;
    effOrd = max(EncOrd.*(EncSnr>=10^(24/20)),[],2) ;
    swtSpt = 1/exp(1) * (2*effOrd+1) ./ (2*pi*snrFrq/340) ;
    trnIdx = find(swtSpt<.0875,1,'first') ;
    trnFrq = snrFrq(trnIdx) ;
    % Generate the decoding option structure
    decOpt = GenerateHoa2SpkOpt('sampFreq',encOpt.sampFreq, ...
        'filterLength',encOpt.filterLength,'transFreq',trnFrq) ;
end

% Decoding filters if necessary
if isempty(decFlt)
    decFlt = Hoa2SpkDecodingFilters(hoaFmtDec,spkFmt,decOpt) ;
    decFlt = decFlt.filters ;
end

% encoding hoaFmt to decoding hoaFmt conversion matrix
Enc2Dec = Hoa2HoaConversionMatrix(hoaFmtEnc,hoaFmtDec) ;


%-----------------------------------%
% FORMATING THE ENC AND DEC FILTERS %
%-----------------------------------%

% Encoding and decoding filters are modified if they don't have the same 
% length and or sampling frequency. First, only the filters with the
% highest sampling frequency are kept.
[smpFrqEnc,smpFrqEncIdx] = max(encOpt.sampFreq) ;
encFlt = encFlt(smpFrqEncIdx) ;
if ~isempty(decOpt.sampFreq)
    [smpFrqDec,smpFrqDecIdx] = max(decOpt.sampFreq) ;
    decFlt = decFlt(smpFrqDecIdx) ;
end

% If encoding and decoding filters have a different sampling frequency, 
% they are resampled to the lowest value. 
if ~isempty(decOpt.sampFreq)
    if smpFrqEnc < smpFrqDec
        firMatrix = zeros(round(smpFrqEnc/smpFrqDec*decOpt.filterLength), ...
            decFlt.nbOutput,decFlt.nbInput) ;
        for I = 1 : hoaFmtDec.nbComp
            firMatrix(:,:,I) = resample(decFlt.firMatrix(:,:,I), ...
                smpFrqEnc,smpFrqDec) ;
        end
        decFlt.firMatrix = firMatrix ;
        decFlt.sampFreq = smpFrqEnc ;
        decOpt.filterLength = size(firMatrix,1) ;
    elseif smpFrqEnc > smpFrqDec
        switch encOpt.filterType
            case 'firMatrix'
                firMatrix = ...
                    zeros(round(smpFrqDec/smpFrqEnc*encOpt.filterLength), ...
                    encFlt.nbOutput,encFlt.nbInput) ;
                for I = 1 : hoaFmtDec.nbComp
                    firMatrix(:,:,I) = resample(encFlt.firMatrix(:,:,I), ...
                        smpFrqDec,smpFrqEnc) ;
                end
                encFlt.firMatrix = firMatrix ;
                encOpt.filterLength = size(firMatrix,1) ;
            case 'postEqFir'
                encFlt.postEqFir = resample(encFlt.postEqFir, ...
                    smpFrqDec,smpFrqEnc) ;
                encOpt.filterLength = size(encFlt.postEqFir,1) ;
        end
        encFlt.sampFreq = smpFrqDec ;
    end
    smpFrq = min(smpFrqEnc,smpFrqDec) ;
else
    smpFrq = smpFrqEnc ;
end

% If some filters are shorter than the others, they are zero-padded.
if ~isempty(decOpt.sampFreq)
    if encOpt.filterLength > decOpt.filterLength
        nmbZerBfr = round((encOpt.filterLength-decOpt.filterLength)/2) ;
        nmbZerAft = encOpt.filterLength-decOpt.filterLength-nmbZerBfr ;
        decFlt.firMatrix = ...
            [ zeros(nmbZerBfr,decFlt.nbOutput,decFlt.nbInput) ; ...
              decFlt.firMatrix ; ...
              zeros(nmbZerAft,decFlt.nbOutput,decFlt.nbInput) ] ;
    elseif encOpt.filterLength < decOpt.filterLength
        nmbZerBfr = round((decOpt.filterLength-encOpt.filterLength)/2) ;
        nmbZerAft = decOpt.filterLength-encOpt.filterLength-nmbZerBfr ;
        switch encOpt.filterType
            case 'firMatrix'
                encFlt.firMatrix = ...
                    [ zeros(nmbZerBfr,encFlt.nbOutput,encFlt.nbInput) ; ...
                      encFlt.firMatrix ; ...
                      zeros(nmbZerAft,encFlt.nbOutput,encFlt.nbInput) ] ;
            case 'postEqFir'
                encFlt.postEqFir = ...
                    [ zeros(nmbZerBfr,size(encFlt.postEqFir,2)) ; ...
                      encFlt.postEqFir ; ...
                      zeros(nmbZerAft,size(encFlt.postEqFir,2)) ] ;
        end
    end
    fltLng = max(encOpt.filterLength,decOpt.filterLength) ;
else
    fltLng = encOpt.filterLength ;
end


%--------------------------------------%
% PLAYBACK FILTERS FREQUENCY RESPONSES %
%--------------------------------------%

% The playback filters are generated in the frequency domain by "array
% multiplicating" the FFTs of encoding and decoding filters.
switch encOpt.filterType
    case 'firMatrix'
        encFft = fft(fftshift(encFlt.firMatrix,1)) ;
    case 'postEqFir'
        encFft = fft(fftshift(encFlt.postEqFir,1)) ;
        encFft = repmat(encFft(:,hoaFmtEnc.index(:,1)+1),[1 1 micFmt.nbMic]) ...
            .* repmat(permute(encFlt.gainMatrix,[3 1 2]),[fltLng 1 1]) ;
end
encFft = permute(encFft,[2 3 1]) ;
plBackFlt = zeros(spkFmt.nbSpk,micFmt.nbMic,fltLng) ;
if ~isempty(decOpt.sampFreq)
    decFft = fft(fftshift(decFlt.firMatrix,1)) ;
    decFft = permute(decFft,[2 3 1]) ;
    for I = 1 : fltLng/2+1
        plBackFlt(:,:,I) = decFft(:,:,I) * Enc2Dec * encFft(:,:,I) ; 
    end
else
    for I = 1 : fltLng/2+1
        plBackFlt(:,:,I) = decFlt.gainMatrix * Enc2Dec * encFft(:,:,I) ; 
    end
end


%------------------------------%
% PLAYBACK FILTER EQUALIZATION %
%------------------------------%

% frequency values
frq = [smpFrq/fltLng/2;(1:fltLng/2)'*smpFrq/fltLng] ;

% A set of source positions
if any(abs(spkFmt.sphCoord(:,2))>1e-3)
    [azmSou,elvSou] = RegularPolyhedron(642) ;
else
    azmSou = (-180:2:178)'*pi/180 ;
    elvSou = zeros(size(azmSou)) ;
end
radSou = Inf(size(azmSou)) ;
souFmt = GenerateSpkFmt('sphCoord',[azmSou elvSou radSou]) ;

% A set of points in the center of the loudspeaker array
radPts = .2 ;
if any(abs(spkFmt.sphCoord(:,2))>1e-3)
    stpPts = .04 ;
    xPts = (-radPts:stpPts:radPts)' ;
    yPts = xPts ;
    zPts = xPts ;
    [xPts,yPts,zPts] = meshgrid(xPts,yPts,zPts) ;
    xyzPts = [ xPts(:) yPts(:) zPts(:) ] ;
    xyzPts = xyzPts(sqrt(sum(xyzPts.^2,2))<=radPts,:) ;
else
    stpPts = .02 ;
    xPts = (-radPts:stpPts:radPts)' ;
    yPts = xPts ;
    [xPts,yPts] = meshgrid(xPts,yPts) ;
    xyzPts = [ xPts(:) yPts(:) zeros(size(xPts(:))) ] ;
    xyzPts = xyzPts(sqrt(sum(xyzPts.^2,2))<=radPts,:) ;
end

% Equalise the filters so that the average energy of the reconstructed 
% sound field is 1 for every frequency
for I = 1 : fltLng/2+1
    
    % The transfers between the sources and the mics
    SouMic = Spk2MicTransfers(souFmt,micFmt,frq(I),1) ;

    % The transfers between the speakers and the points
    SpkPts = SphericalWave(spkFmt.xyzCoord,xyzPts,frq(I),1) ;
    
    % Reconstructed sound field at the points for every source direction
    SndFld = SpkPts * plBackFlt(:,:,I) * SouMic ;
        
    % Playback filter equalization
    plBackFlt(:,:,I) = plBackFlt(:,:,I) / norm(SndFld,'fro') ;
    
end


%------------------------------------%
% PLAYBACK FILTERS IMPULSE RESPONSES %
%------------------------------------%

% Back in the time domain
plBackFlt = permute(plBackFlt,[3 1 2]) ;
plBackFlt(fltLng/2+2:end,:,:) = conj(plBackFlt(fltLng/2:-1:2,:,:)) ;
plBackFlt = fftshift(real(ifft(plBackFlt)),1) ;

% Time-window the filters
plBackFlt = repmat(hanning(fltLng),[1 spkFmt.nbSpk micFmt.nbMic]) ...
    .* plBackFlt ;

% Normalise the filters so that the maximum amplification is 0dB
plBackFlt = plBackFlt / max(max(max(abs(fft(plBackFlt))))) ;




