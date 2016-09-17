function [preRep,preOrg,xPts,yPts] = ...
    Mic2Hoa2SpkPressureAroundHead(mic2HoaCfg,D,spkSphCoord,freq,plotMode,souSphCoord)

% Plots the orignal pressure field obtained around a rigid sphere with a
% radius of 10cm when a single source located at souSphCoord plays at 
% frequency freq.
% Compare with the sound field obtained when:
% - recording/encoding the same source using microphone/encoding filters
%   defined in mic2HoaCfg
% - decoding the HOA signals using D decoding matrix and playing back 
%   over the loudspeakers located at spkSphCoord

% Default value for the source spherical coordinates
if nargin < 6
    souSphCoord = [0 0 2] ;
end

% Default visualisation mode ('real' or 'comp')
if nargin < 5
    plotMode = 'comp' ;
end

% Checks if the decoding matrix is coherent with the speaker array
if size(D,1) ~= size(spkSphCoord,1)
    error('The decoding matrix is not coherent with the loudspeaker coordinates') ;
end

% Checks if the decoding matrix is coherent with the encoding structure
if size(D,2) ~= mic2HoaCfg.hoaFmt.nbComp
    error('The decoding matrix is not coherent with the encoding structure') ;
end

% To make the code easier, structures are separated
hoaFmt     = mic2HoaCfg.hoaFmt ;
micFmt     = mic2HoaCfg.micFmt ;
mic2HoaOpt = mic2HoaCfg.mic2HoaOpt ;
filters    = mic2HoaCfg.filters ;

% Encoding error is only computed for the highest sampling frequency
% value available
[fs,indFs] = max(mic2HoaOpt.sampFreq) ;
fftNb = mic2HoaOpt.filterLength(indFs) ;
filters = filters(indFs) ;

% Just a little test
if max(freq) > fs/2
    error('f cannot be higher than half the sampling frequency') ;
end


%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALISATION POINTS %
%%%%%%%%%%%%%%%%%%%%%%%%

% azimuth and radius resolution (in deg and mm)
resAzm = 1 ;
resRad = 2.5 ;

% azimuth and radius values
azmPts = (0:resAzm:360)'*pi/180 ;
radPts = (.1:resRad/1000:.15)' ; 
[azmPts,radPts] = meshgrid(azmPts,radPts) ;

% cartesian coordinates
[xPts,yPts] = sph2cart(azmPts,zeros(size(azmPts)),radPts) ;

% "MicFmt" corresponding to the visualisation points
micFmtPts = GenerateMicFmt( ...
    {'sphCoord',[azmPts(:),zeros(size(azmPts(:))),radPts(:)]},.1) ;

% HOA decomposition used for modelling the problem
hoaOrdMod = Hoa2MicRequiredOrder(micFmtPts,freq,max(hoaFmt.index(:,1))) ;
if ~(any(spkSphCoord(:,2))~=0)&&(souSphCoord(:,2)==0)
    hoaFmtMod = GenerateHoaFmt('res2d',hoaOrdMod,'res3d',hoaOrdMod, ...
        'plane','xOy') ;
else
    hoaFmtMod = GenerateHoaFmt('res2d',hoaOrdMod,'res3d',hoaOrdMod) ;
end

% Hoa -> Pts transfer matrix 
Tpts = Hoa2MicTransferMatrix(hoaFmtMod,micFmtPts,freq) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORIGINAL SOUNDFIELD HOA DECOMPOSITION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Original pressure = 1 at the center
Wsou = Spk2HoaRadialWeightings(0:hoaOrdMod, ...
    2*pi*freq/340*souSphCoord(3)) ;
Ysou = SphericalHarmonicMatrix(hoaFmtMod, ...
    souSphCoord(:,1),souSphCoord(:,2)) ;
b = Wsou(hoaFmtMod.index(:,1)+1) .* Ysou ;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% (HOA)->(MIC) TRANSFERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%

Tmic = Hoa2MicTransferMatrix(hoaFmtMod,micFmt,freq) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (MIC)->(ESTIMATED-HOA) TRANSFERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Encoding filter ferquency responses
switch mic2HoaOpt.filterType
    case 'postEqFir'
        % Equalization filter Fourier Transforms are computed
        postEqFir = fft(fftshift(filters.postEqFir,1),[],1) ;
        postEqFir = postEqFir(:,hoaFmt.index(:,1)+1) ;
        % Global encoding filter responses (including gain matrix)
        FltRsp = repmat(postEqFir,[1 1 micFmt.nbMic]) ...
            .* repmat(permute(filters.gainMatrix,[3 1 2]),[fftNb 1 1]) ;
    case 'firMatrix'
        % Encoding filter Fourier Transforms are computed
        FltRsp = fft(fftshift(filters.firMatrix,1),[],1) ;
end

% Encoding filter responses at 'freq' 
fftFreq = (0:fftNb-1).'*fs/fftNb ;
FltRsp = interp1(fftFreq,FltRsp,freq) ;
FltRsp = permute(FltRsp,[2 3 1]) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% (SPK)->(HOA) TRANSFERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% speaker radial weightings
Wspk = Spk2HoaRadialWeightings(0:hoaOrdMod, ...
    2*pi*freq/340*spkSphCoord(:,3)) ;

% speaker spherical harmonic values
Yspk = SphericalHarmonicMatrix(hoaFmtMod, ...
    spkSphCoord(:,1),spkSphCoord(:,2)) ;

% Spk2Hoa transfers
Tspk = Wspk(hoaFmtMod.index(:,1)+1,:) .* Yspk ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORIGINAL AND REPRODUCED SOUND FIELDS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Original sound field
preOrg = reshape(Tpts*b,size(xPts)) ;

% Reproduced soundfield
preRep = reshape(Tpts*Tspk*D*FltRsp*Tmic*b,size(xPts)) ;

%%%%%%%%%
% PLOTS %
%%%%%%%%%

switch plotMode
    
    case 'real'
        
        % Reproduced sound field
        figure ;
        surf(xPts,yPts,real(preRep)) ;
        view(2) ; axis equal ; axis off ;
        colormap gray ; caxis([-1 1]) ; colorbar ; shading interp ;
        title('Reproduced')
        
        % Original sound field
        figure ;
        surf(xPts,yPts,real(preOrg)) ;
        view(2) ; axis equal ; axis off ;
        colormap gray ; caxis([-1 1]) ; colorbar ; shading interp ;
        title('Original')

    case 'comp'
        
        % Sound field amplitude
        figure ;
        surf(xPts,yPts,20*log10(abs(preRep))) ;
        view(2) ; axis equal ; axis off ;
        caxis([-12 12]) ; colorbar ; shading interp ;
        title('Reproduced: amplitude')
        figure ;
        surf(xPts,yPts,20*log10(abs(preOrg))) ;
        view(2) ; axis equal ; axis off ;
        caxis([-12 12]) ; colorbar ; shading interp ;
        title('Original: amplitude')
        
        % Sound field phase
        figure ;
        surf(xPts,yPts,angle(preRep)) ;
        view(2) ; axis equal ; axis off ;
        colormap hsv ; caxis([-pi pi]) ; colorbar ; shading interp ;
        title('Reproduced: phase')
        figure ;
        surf(xPts,yPts,angle(preOrg)) ;
        view(2) ; axis equal ; axis off ;
        colormap hsv ; caxis([-pi pi]) ; colorbar ; shading interp ;
        title('Original: phase')

end



