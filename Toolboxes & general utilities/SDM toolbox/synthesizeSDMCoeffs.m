function [H,Hbin] = synthesizeSDMCoeffs(P,DOA, p)
% H = synthesizeSDMCoeffs(P,DOA,varargin)
% Synthesizes the pressure vector P according to DOA matrix 
% to loudspeakers defined in p.lspLocs
% 
% USAGE:
% 
% H : Synthesized matrix of impulse responses [N numberOfLoudspeakers]
% Hbin : Synthesized binaural impulse responses [2 numberOfLoudspeakers]
%        Channels 1: left, 2: right
% 
% P : Vector with the pressure in the center of the array [N 1]   
% DOA : Matrix DOA corresponding to each pressure value [N 3]
% p: struct given by createSynthesisStruct.m or can be manually given
% 
% EXAMPLES:
% 
% For further examples for setting the parameters, see 
% createSynthesiStruct.m or demo*.m
% 
% References
% [1] S. Tervo, J. Pätynen, A. Kuusinen, T. Lokki "Spatial decomposition method 
% for room impulse responses", Journal of the Audio
% [2] S. Tervo, J. Pätynen, N. Kaplanis, M. Lydolf, S. Bech, and T. Lokki
% "Spatial Analysis and Synthesis of Car Audio System and Car Cabin
% Acoustics with a Compact Microphone Array",

% SDM toolbox : synthesizeSDMCoeffs
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Sakari.Tervo@aalto.fi and Jukka.Patynen@aalto.fi

if nargin < 3
    error(['SDM toolbox : synthesizesSDMCoeffs : Three inputs are required ' ... 
    ' synthesizesSDMCoeffs(P,DOA, p)']);
end

disp('Started synthesizing'); ts = tic;

% Local variable for loudspeaker locations
lspLocs = p.lspLocs;

numOfLsp = size(lspLocs,1);

% Assign NaNs to the LEF channels
lspLocs(p.LFEchannel,:) = NaN;

% To spherical coordinates
[az,el] = cart2sph(DOA(:,1), DOA(:,2), DOA(:,3));

% ---- hack ----
% Sometimes you get NaNs from the DOA analysis
% Replace NaN directions with uniformly distributed random angle
az(isnan(az)) = pi-rand(size(az(isnan(az))))*2*pi;
el(isnan(el)) = pi/2-rand(size(el(isnan(el))))*pi;
% ------ EOF hack -------

% Nearest loudspeaker synthesis == no panning
[nLOC(:,1),nLOC(:,2),nLOC(:,3)] = sph2cart(az,el,1);
[lpsx,lpsy,lpsz] = sph2cart(lspLocs(:,1)/180*pi, lspLocs(:,2)/180*pi, 1);

% Nearest neighbour search
indslps = knnsearch([lpsx,lpsy,lpsz],nLOC);

% Output matrix of impulse responses for each loudspeaker
H = zeros(length(P),numOfLsp);

% p.dimensionality is user defined, typically 2-D arrays are in home use,
% and 3-D arrays in research facilities and cinemas (movie theaters)

% ----- Synthesis of the impulse response as in Ref. [1] ----------------
% ---- For a 3-D loudspeaker setup, e.g. Hamasaki 22.2 -----
if p.dimensionality == 3  
    % assign the pressure sample to nearest loudspeakers w.r.t. az and el
    for rep_lps = 1:numOfLsp
        inds = indslps == rep_lps;
        H(inds,rep_lps) = P(inds); % assign the pressure values
    end
    
elseif p.dimensionality == 2
    % ---- For a 2-D loudspeaker setup (azimuthal plane) -----
    % The pressure is weighted with the elevation, thus all the energy is played back
    % but higher elevation is attenuated more (according to abs(cos(el))).
    
    % Assign the pressure sample to nearest loudspeakers w.r.t. az and el
    for rep_lps= 1:numOfLsp
        inds = indslps == rep_lps;
        H(inds,rep_lps) = P(inds).*abs(cos(el(inds)));
    end
end
% ----> H is the synthesized impulse response

% The rest is some tricks to obtain a smoother time-frequency response

% ------- Smoothen the spectrogram of H to match P ----------
% Post-equalization of the synthesized impulse responses
% This method is described in Reference [2]
% FIR octave band filters
HS = zeros(size(H));
snfft = length(P);
numOfBands = length(p.f1);
% ---- post-eq octave bands -----
for band = 1:numOfBands
    if band == 1
        winLen = snfft;
        winLen = winLen + mod(winLen,2);
    else
        winLen = round(7/p.f1(band)*p.fs);
        winLen = winLen + mod(winLen,2);      
    end
    % Equalize the NSL output in winLen sized frames
    H_posteq = equalizeNLS(H, P, winLen);
    
    % Filter the result with octave band filter
    H_filt = fftconv(p.g(:,band), H_posteq, snfft);
    HS = HS + H_filt;
end
H = HS; % overwrite H with the post-equalized version


% Since there are several loudspeakers, the low-frequencies
% are low-passed due to phase summation in the listening location.
% Here, the cut-off frequency is 200 Hz, which assumes that in-phase
% summation occurs below 200 Hz
if numel(p.LFEchannel) == 0
    f_cutoff = 200;
    [b,a] = butter(1,f_cutoff/(p.fs),'high');
    H = filtfilt(b,a,H);
end

% ----------- Add LFE channels --------------
% If one LFE channel
if numel(p.LFEchannel) == 1
    
    % High-passed sound to all other channels
    Hhp = filtfilt(p.Bhp,p.Ahp,H);
    H = Hhp;
    
    % Low-passed sound to LFE channel
    Hlp = filtfilt(p.Blp,p.Alp,H);
    H(:,p.LFEchannel) = sum(Hlp,2);
    
% If two LFE channels
elseif numel(p.LFEchannel) == 2 % This assumes that LFE1 is left and LFE2 is right 
    
    % High-passed sound to all other channels
    Hhp = filtfilt(p.Bhp,p.Ahp,H);
    H = Hhp;
    
    % Low-passed sound to LFE channels
    inds_lfe1 = lspLocs(:,1) > 0;
    inds_lfe2 = lspLocs(:,1) < 0;
    inds_zero = lspLocs(:,1) == 0;
    
    Hlp1 = filtfilt(p.Blp,p.Alp,[H(:,inds_lfe1), sqrt(2)*H(:,inds_zero)]);
    Hlp2 = filtfilt(p.Blp,p.Alp,[H(:,inds_lfe2), sqrt(2)*H(:,inds_zero)]);

    H(:,p.LFEchannel(1)) = sum(Hlp1,2);
    H(:,p.LFEchannel(2)) = sum(Hlp2,2);
end
% TODO, elseif more than two LFE channels, depends on the positions of the other
% loudspeakers

%  ------ Binaural synthesis  ---------------
if p.Binaural
    disp('Started binaural synthesis.');tb = tic;
    N1 = size(H,1);
    N2 = size(p.hrir_l,1);
    nfft = N1+N2+1;
    Hbin(:,1) = real(ifft(sum(fft(p.hrir_l,nfft).*fft(H,nfft),2)));
    Hbin(:,2) = real(ifft(sum(fft(p.hrir_r,nfft).*fft(H,nfft),2)));
    disp(['Ended binaural synthesis in ' num2str(toc(tb)) ' seconds.' ])
end

% --- Compensate for the difference between the distances of the
% loudspeakers ----
delays = round((max(p.lspLocs(:,3)) - p.lspLocs(:,3))/p.c*p.fs);

for lsp = 1:numOfLsp
    H(:,lsp) = [zeros(delays(lsp),1); H(1:end-delays(lsp),lsp)];   
end

% --- Draw the array if asked -----
if p.showArray
    plotArray(p.lspLocs, p.LFEchannel)
end

disp(['Ended synthesizing in ' num2str(toc(ts)) ' seconds.' ])

end

function plotArray(lspLocs, LFEchannel)
% plotArray(micLocs)
% Plot the geometry of loudspeaker locations, defined in lspLocs

% SDM toolbox : plotArray
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Copyleft

[lspLocs(:,1),lspLocs(:,2), lspLocs(:,3)] = ...
    sph2cart(lspLocs(:,1)/180*pi,lspLocs(:,2)/180*pi, lspLocs(:,3));

numOfLsp = size(lspLocs,1);
figure
hold on
lfe = 1;
for m = 1:numOfLsp
    if any(m == LFEchannel)
        plot3(lspLocs(m,1),lspLocs(m,2),lspLocs(m,3),'color','k','Marker','s')
        text(lspLocs(m,1),lspLocs(m,2),lspLocs(m,3)-0.05,[ '\# ' num2str(m) ' LFE' num2str(lfe)]);
        lfe = lfe + 1;
    else
        plot3(lspLocs(m,1),lspLocs(m,2),lspLocs(m,3),'color','k','Marker','o')
        text(lspLocs(m,1),lspLocs(m,2),lspLocs(m,3),[ '\# ' num2str(m)]);
        
    end
end

hold off;grid on;
view([35 45]);
axis square;
title('Loudspeaker Array Geometry in the Synthesis');
xlabel('X-coord.');ylabel('Y-coord.');zlabel('Z-coord.');

end

function xm_comp = equalizeNLS(x_sdm,x_ref,winLen)
% xm_comp = equalizeNLS(x_sdm,x_ref,winLen)
% Matches the x_sdm spectrogram to x_ref with overlap-add in a windows of
% size winLen.
%
% USAGE:
% 
% x_comp : Spectrogram matched impulse responses, same size as x_sdm
% x_sdm  : Matrix of impulse responses subject to compensation [N numOfLoudspeakers]
% x_ref  : Vector including reference impulse response for compensation [N 1]
% winLen : window length in samples
% 
% EXAMPLES:
% 
% xm_comp = equalizeNLS(randn(1000,24),randn(1000,1),128);
%
% For further usage see See synthesizesSDMCoeffs.m

% SDM toolbox : equalizeNLS
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Copyleft

% Add zeros to fully construct the IR
N = size(x_sdm,1);
x_sdm = [zeros(winLen,size(x_sdm,2));x_sdm;zeros(winLen,size(x_sdm,2))];
x_ref = [zeros(winLen,size(x_ref,2));x_ref;zeros(winLen,size(x_ref,2))];

% Overlap-add parameters
win = hanning(winLen);
win = win/max(win);
nfft = winLen;
startx = 1;
endx = winLen;
nforward = winLen/2;

% Number of frames, i.e., windows
NN = round((size(x_sdm,1)-winLen)/nforward);

y = zeros(size(x_sdm,1)+2*nfft, size(x_sdm,2));
for n = 1:NN
    
    % Disply the current frame number
    if(mod(n,1000)==0)
        disp(['equalizeNLS: processing frame : ' num2str(n)])
    end
    
    % Reference frame and NLS impulse response frames
    X_ref = fft(bsxfun(@times, win, x_ref(startx:endx)),nfft*2);
    X_sdm = fft(bsxfun(@times, win, x_sdm(startx:endx,:)),nfft*2);
    
    % Compute the difference
    comp = sqrt(abs(X_ref).^2 ./ (sum(abs(X_sdm).^2,2)+eps));
    % Square
    Hf = sqrt(comp).*sqrt(conj(comp));
    % Calculate the compensated version of the NLS frames
    Yy = bsxfun(@times, Hf, X_sdm);
    % Back to the time-domain
    tempy = real(ifft(Yy));
    
    % Overlap-add
    y(startx:startx + nfft*2-1,:) = y(startx:startx + nfft*2-1,:) +  tempy;
    startx = startx + nforward;
    endx = endx + nforward;
    
end
% Truncate to the original size
xm_comp = y(winLen+1:winLen+N,:);

end

function y = fftconv(b, x, snfft)
% Convolution via fft for two impulse responses b and x
% b : First FIR
% x : Another FIR
% snfft : length of the result

% SDM toolbox : equalizeNLS
% Sakari Tervo & Jukka Pätynen, Aalto University, 2011-2016
% Copyleft

y = real(ifft(bsxfun(@times, fft(b,2*snfft), fft(x,2*snfft))));

y = y(snfft+1:end,:);
y = y(1:snfft,:);
end