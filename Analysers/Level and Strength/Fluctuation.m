function [OUT, varargout] = Fluctuation(IN,soundfield,signaltype)

% This function attempts to implement a novel fluctuation strength model 
% for modulated and unmodulated signals in free and diffuse sound fields as 
% detailed in Tingting Zhou et al. "A Model for Calculating 
% Psychoacoustical Fluctuation Strength", JAES, September 2015.
%
% The model uses 75 ERB channels to calculate generalized modulation depth 
% (GMD) values, which are converted into specific fluctuation strengths.
% Test signal: 60dB 1kHz 100% mod 4Hz = 1 vacil
%
% INPUT ARGUMENTS:
% IN - audio signal (1 channel only)
% soundfield - free field = 0; diffuse field = 1
% signaltype - unmodulated = 0; modulated = 1;
%
% OUTPUTS:
% Total Fluctuation strength
% Fluctuation strength statistics
% Time-averaged specific Fluctuation strength
% Time-varying Fluctuation strength
%
% MATLAB code by Ella Manor for AARAE (2015) 

% *************************************************************************

if isstruct(IN)
        
    IN = choose_from_higher_dimensions(IN,1,1); % only 1 channel analysis at present
    audio = IN.audio;
    fs = IN.fs;      
    
    if isfield(IN,'cal')
        cal = IN.cal;
    else
        warndlg('Calibration data missing - please calibrate prior to calling this function.','AARAE info','modal');
        IN = cal_aarae(IN);
        cal = IN.cal;
    end

    if isfield(IN,'name') 
        name = IN.name;
    else
        name = '';
    end
else
    audio = IN(:,1,1,1,1,1);
    cal = cal(1);
    name = '';
end

% *************************************************************************

if nargin == 1 

    param = inputdlg({'Sound field: Free (0); Diffuse (1)';... 
                      'Sound signal: Unmodulated (0); Modulated (1)'},...
                      'Parameters',... 
                      [1 50],... 
                      {'0';'1'});
    
    if length(param) < 2, param = []; end 
                                          
                                         
    if ~isempty(param) 
                      
        soundfield = str2double(param(1));
        signaltype = str2double(param(2));
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
    
else
    param = [];
end

% *********************************************************************

if ~isempty(audio) && ~isempty(fs)
    
    calconstant = 0;
    cal = cal-calconstant;
    if length(cal) == 1
        audio = audio .* 10.^(cal/20);
    else
        audio(:,1) = audio(:,1) .* 10.^(cal(1)/20);
        audio(:,2) = audio(:,2) .* 10.^(cal(2)/20);
    end

    disp(['rms level of the entire wave ', num2str(10*log10(mean(audio.^2)+10e-99)+calconstant), ' dB'])

[samples,nchan] =  size(audio);

    if nchan > 1
        audio = audio(:,1);
        disp('Only one channel signals are analysed by Loudness_ISO532.') 
        nchan = 1;
    end
    
% *********************************************************************
N = 1024;
window = blackman(N);

% ERB filter design is adopted from MGBLoudness2b.m
% roex
f = (1:N)*fs/N;
nfreq = length(f);
erb = 24.673.*(0.004368.*f+1);
p = 4.*f./erb;
g = abs((repmat(f',1,nfreq)-repmat(f,nfreq,1))./repmat(f,nfreq,1)); 
[row,col] = find(g <= 2);
indg = find(g <= 2);

% excitation
p511k = 30.2012922; %p51 @ 1kHz - see Glasberg and Moore (1990)
ERBS = 1.8:0.5:38.9;
nERBS = length(ERBS);
ERBSfreq = 1000*(exp((ERBS./21.36554).*log(10))-1)./4.368;
ERBSerb = 24.673 * (4.368 .* ERBSfreq./1000 + 1);
p51 = 4 .* ERBSfreq./ERBSerb;
g2 = (repmat(f',1,nERBS)-repmat(ERBSfreq,nfreq,1))./repmat(ERBSfreq,nfreq,1);
p2 = zeros(nfreq,nERBS);
g2lessthan2 = g2<2;
g2pos = (g2>=0) & g2lessthan2;
p51=repmat(p51,nfreq,1);
p2(g2pos) = p51(g2pos);
g2neg = g2<0;
[row_g2_neg, ~]=find(g2<0);

% pad audio signal for window base processing
if mod(samples,N) ~= 0
    n = ceil(samples/N)+1;
else
    n = ceil(samples/N);
end
nn = n*N-samples;
if mod(nn,2)
    padsamp = zeros(nn/2,nchan);
else
    padsamp = zeros(nn/2,nchan);
end
audio = [padsamp;audio;padsamp];
samples = length(audio);

% hearing threshold per Bark number (adopted from RoughnessDW.m)
HTres= [0 130
        0.01 70
        0.17 60
        0.8 30
        1 25
        1.5	20
        2 15
        3.3	10
        4 8.1
        5 6.3
        6 5
        8 3.5
        10 2.5
        12 1.7
        13.3 0
        15 -2.5
        16 -4
        17 -3.7
        18 -1.5
        19 1.4
        20 3.8
        21 5
        22 7.5
        23 15
        24 48
        24.5 60
        25 130];

% hearing threshold values interpolated for window length

NN=linspace(0,25,N);
MinExcdB = interp1(HTres(:,2),NN','linear','extrap');

% Transfer functions for the outer and middle ear.
if soundfield == 0 % free field
% Free field to eardrum transfer function (frontal incinece). First column
% = frequnecy band in Hz. Second column = level at the eardrum - level
% measured in free field in dB [ANSI S3.4-2007]
filt_outer = [0 0.;
                20 0.;
                25 0.;
                31.5 0.;
                40 0.;
                50 0.;
                63 0.;
                80 0.;
                100 0.;
                125 0.1;
                160 0.3;
                200 0.5;
                250 0.9;
                315 1.4;
                400 1.6;
                500 1.7;
                630 2.5;
                750 2.7;
                800 2.6;
                1000 2.6;
                1250 3.2;
                1500 5.2;
                1600 6.6;
                2000 12.;
                2500 16.8;
                3000 15.3;
                3150 15.2;
                4000 14.2;
                5000 10.7;
                6000 7.1;
                6300 6.4;
                8000 1.8;
                9000 -0.9;
                10000 -1.6;
                11200 1.9;
                12500 4.9;
                14000 2.0;
                15000 -2.0;
                16000 2.5;
                20000 2.5];

b_outer = fir2(N-1,[filt_outer(:,1)./(0.5*fs);1],[10.^(filt_outer(:,2)./20);0]);

elseif soundfield == 1 % diffuse field

% Diffuse field to eardrum transfer function. First column
% = frequnecy band in Hz. Second column = level at the eardrum - level
% measured in diffuse field in dB [ANSI S3.4-2007]
filt_outer = [0 0.;
                20 0.;
                25 0.;
                31.5 0.;
                40 0.;
                50 0.;
                63 0.;
                80 0.;
                100 0.;
                125 0.1;
                160 0.3;
                200 0.4;
                250 0.5;
                315 1.;
                400 1.6;
                500 1.7;
                630 2.2;
                750 2.7;
                800 2.9;
                1000 3.8;
                1250 5.3;
                1500 6.8;
                1600 7.2;
                2000 10.2;
                2500 14.9;
                3000 14.5;
                3150 14.4;
                4000 12.7;
                5000 10.8;
                6000 8.9;
                6300 8.7;
                8000 8.5;
                9000 6.2;
                10000 5.;
                11200 4.5;
                12500 4.0;
                14000 3.3;
                15000 2.6;
                16000 2.;
                20000 2.];

b_outer = fir2(N-1,[filt_outer(:,1)./(0.5*fs);1],[10.^(filt_outer(:,2)./20);0]); 

end
% Transfer function of middle ear. First column
% = frequnecy band in Hz. Second column = effective level at the cochlea
% relative to the level at the eardrum in dB. Values of 14kHz-20kHz have
% not been validated. 
filt_middle = [0 0.;
            20 -39.6;
            25 -32.;
            31.5 -25.85;
            40 -21.4;
            50 -18.5;
            63 -15.9;
            80 -14.1;
            100 -12.4;
            125 -11.;
            160 -9.6;
            200 -8.3;
            250 -7.4;
            315 -6.2;
            400 -4.8;
            500 -3.8;
            630 -3.3;
            750 -2.9;
            800 -2.6;
            1000 -2.6;
            1250 -4.5;
            1500 -5.4;
            1600 -6.1;
            2000 -8.5;
            2500 -10.4;
            3000 -7.3;
            3150 -7.;
            4000 -6.6;
            5000 -7.;
            6000 -9.2;
            6300 -10.2;
            8000 -12.2;
            9000 -10.8;
            10000 -10.1;
            11200 -12.7;
            12500 -15.;
            14000 -18.2;
            15000 -23.8;
            16000 -32.3;
            18000 -45.5
            20000 50.];            

b_middle = fir2(N-1,[filt_middle(:,1)./(0.5*fs);1],[10.^(filt_middle(:,2)./20);0]);

% Filtering of the audio signal in the frequency domain with transfer 
% functions for outer and middle ear. Spectrum components bellow the 
% hearing threshold are ignored. 
startIndex = 1;
endIndex = N;
filteredaudio = zeros(samples,nchan);
SPL_mat = zeros(1,n); % sound pressure level of signal before processing
TimePoints = zeros(n,1);                                    
for windowNum = 1:n
    dataIn = audio(startIndex:endIndex,nchan);
    
    SPL = mean(rms(dataIn));
        if SPL > 0
            SPL = mag2db(SPL);
        else
            SPL = -400;
        end
    SPL_mat(windowNum) = SPL;

    SpectrumIn = fft(dataIn.*window);
    OuterIn = filter(b_outer,1, SpectrumIn);
    MiddleIn = filter(b_middle,1, OuterIn);
    Lg = abs(MiddleIn);
    LdB	= mag2db(Lg);
    whichL = find(LdB>MinExcdB);
    sizL = length(whichL);
    etmp = zeros(1,N);
    for l = 1:1:sizL
        etmp(whichL(l)) = MiddleIn(whichL(l),nchan);    
    end   
    filteredaudio(startIndex:endIndex,nchan) = etmp;
    startIndex = startIndex+N;
    endIndex = endIndex+N;
    TimePoints(windowNum,1) = startIndex/fs; 
end

% Conver signal from frequency to time domain
filteredTemp = abs(N*real(ifft(filteredaudio)));

% Calculate excitation levels for each of the 75 ERB channels, adopted from
% MGBLoudness2b.m
startIndex = 1;
endIndex = N;
Elevel = zeros(nERBS,samples);
Exsig = zeros(nERBS,samples);
for windowNum = 1:n
    dataIn = filteredTemp(startIndex:endIndex,nchan);
    intensity = zeros(nfreq,nfreq);
    intensity(indg) = (1+p(col)'.*(g(indg))).* exp(-p(col)'.*(g(indg))) .* dataIn(row,1);
    RoexL = 10.*log10(sum(intensity,1)+eps);
    
    p2(g2neg) = p51(g2neg)-0.35.*(p51(g2neg)./p511k) .* (RoexL(row_g2_neg(:))-51)';
    if p2(g2lessthan2) <0.1, p2(g2lessthan2) = 0.1; end
    
    Imat = repmat(dataIn,1,nERBS);
    Exsig(:,startIndex:endIndex) = (g2lessthan2.*Imat .* (1+p2.*abs(g2)).*exp(-p2.*abs(g2)))'+eps;
    Elevel(:,startIndex:endIndex) = 10.*log10(Exsig(:,startIndex:endIndex));
    startIndex = startIndex+N;
    endIndex = endIndex+N;
end

% Calculate generalised modulation depth (GMD) for each channel
% for modulated signal:
if signaltype == 1 
    transform = hilbert(Exsig);
%     transformreal = real(transform);
    transformimag = imag(transform); % the transform
    env = sqrt(Exsig.^2 + transformimag.^2);
    
% for unmodulated signals the envelope is extracted by finding the local
% maxima values
elseif signaltype == 0
    env = max(Exsig,1);
end

% Bandpass weighting function 0.25 - 64 Hz as described in the paper, Fig.1
cf1 = 0.25./(0.5*fs);
cf2 = 64./(0.5*fs);
b_filt_bandpass = fir1(N,[cf1 cf2],'bandpass');

% Filter the signal with the filter
startIndex = 1;
endIndex = N;
for windowNum = 1:n
    dataIn = fft(env(:,startIndex:endIndex));
    hBPi = filter(b_filt_bandpass,1,dataIn);
    hBPit(:,startIndex:endIndex) = abs(N*real(ifft(hBPi)));
    startIndex = startIndex+N;
    endIndex = endIndex+N;
end

% Relative fluctuation strength as a function of Modulation depth (dB),
% Fig. 2
fig2 = [0 0.;
        1. 0.05;
        2 0.12;
        3 0.14;
        4 0.18;
        5 0.2;
        6 0.29;
        7 0.37;
        8 0.42;
        9 0.51;
        10 0.6;
        20 0.9;
        30 0.98;
        40 1.];

% Weightning coefficient in each ERB channel 
gerb = [0.85 0.91 0.95 1.03 1.1 1.16 1.18 1.2 1.223 1.228,...
    1.233 1.238 1.243 1.248 1.25 1.26 1.27 1.275 1.3 1.32,...
    1.34 1.33 1.34 1.343 1.348 1.35 1.34 1.343 1.35 1.345,...
    1.344 1.343 1.342 1.34 1.338 1.333 1.33 1.32 1.31 1.3,...
    1.29 1.285 1.28 1.27 1.26 1.25 1.235 1.23 1.225 1.21,...
    1.205 1.2 1.17 1.16 1.14 1.11 1.09 1.07 1.05 1.04,...
    1.035 1. 0.99 0.95 0.93 0.92 0.91 0.86 0.84 0.82 0.79 0.76 0.75 0.74 0.71];    

gerbierb = 1:1:75;
gerbi = (interp1(gerb',gerbierb'));

% Calculating specific fluctuation strength for each of the 75 channels
startIndex = 1;
endIndex = N;
mi_mat = zeros(nERBS,n);
fi_mat = zeros(nERBS,n);
wm = zeros(nERBS,n);
% env CONVERT TO TIME DOMAIN
for windowNum = 1:n
    dataInT = env(:,startIndex:endIndex); % envelope before bandpass filtering
    dataIn = hBPit(:,startIndex:endIndex); % envelope after bandpass filtering
    h0 = mean(dataInT,2);
    hBPrms = rms(dataIn,2);
    mi = zeros(1,nERBS);
    mi_depth = zeros(1,nERBS);
    for k = 1:75
        if h0(k)>0
            mi(k) = hBPrms(k)/h0(k);
            if mi(k)>1
                mi(k)=0.99;
            end
        else
            mi(k)=0;
        end
        mi_depth(k) = 20.*log10((1+mi(k)./(1-mi(k))));
        xx = mi_depth(k)-fig2(:,1);
        [~,xi] = min(abs(xx));
        wm(k,windowNum) = fig2(xi,2);
        fi_mat(k,windowNum) = gerb(k).*wm(k,windowNum);
    end
    
    % Calculating cross correlation coefficients between excitation signal
    % envelope of different channels, code adopted from RoughnessDW.m
    ki = zeros(1,nERBS);
    for k = 1:1:73
        cfac = cov(dataIn(k,:),dataIn(k+2,:));
        den	= diag(cfac);
        den	= sqrt(den*den');
        if den(2,1)>0
            ki(k) =	cfac(2,1)/den(2,1);
        else
            ki(k) =	0;
        end
    end
    fi_mat(1,windowNum) = fi_mat(1,windowNum)*(ki(1)^2);
    fi_mat(2,windowNum) = fi_mat(2,windowNum)*(ki(2)^2);
    for k = 3:1:73
        fi_mat(k,windowNum)	= fi_mat(k,windowNum)*((ki(k-2)*ki(k))^2);
    end
    fi_mat(74,windowNum) = fi_mat(74,windowNum)*(ki(72)^2);
    fi_mat(75,windowNum) = fi_mat(75,windowNum)*(ki(73)^2);
    
    mi_mat(:,windowNum) = mi; 
    startIndex = startIndex+N;
    endIndex = endIndex+N;
end

Fi = mean(fi_mat,2); % Time-averaged specific fluctuation strength 
% in each channel
FiT = mean(fi_mat,1); % Time-varying fluctuation strength
Cal = 0.03; % Adjust total fluctuation strength, so that 
% test signal Ftotal ~ 1.0005 vacil
F = Cal*sum(Fi); % Total fluctuation strength

% *************************************************************************
% Data Presentation
% *************************************************************************

if isstruct(IN)
    
    % ********* TABLES *********
    
    % Fluctuation statistics
    Fmean = mean(Fi);
    Fstd = std(Fi);
    Fmax = max(Fi);
    F1 = prctile(Fi,99);
    F2 = prctile(Fi,98);
    F3 = prctile(Fi,97);
    F4 = prctile(Fi,96);
    F5 = prctile(Fi,95);
    F10 = prctile(Fi,90);
    F20 = prctile(Fi,80);
    F30 = prctile(Fi,70);
    F40 = prctile(Fi,60);
    F50 = median(Fi);
    F60 = prctile(Fi,40);
    F70 = prctile(Fi,30);
    F80 = prctile(Fi,20);
    F90 = prctile(Fi,10);
    Fmin = min(Fi);
    
    dataF = [F;Fmean;Fstd;Fmax;F1;F2;F3;F4;F5;F10;F20;F30;F40;F50;F60;F70;F80;F90;Fmin];
    
    % generate tables of results
    
    fig1 = figure('Name','Fluctuation Strength Statistics');
    table1 = uitable('Data',dataF,...
        'ColumnName',{'Fluctuation Strength'},...
        'RowName',{'Total fluctuation','Mean','Standard deviation','Maximum',...
        'F1','F2','F3','F4',...
        'F5','F10','F20','F30','F40','F50 (median)','F60',...
        'F70','F80','F90','Minimum'});
    
    [~,tables] = disptables(fig1,table1); % AARAE function
    
    OUT.tables = tables;
%  
    % ********* CHARTS *********
    
    figure('Name',['Fluctuation Strength ',name])
    
    subplot(2,2,1)
    
    % Time-varying SPL to monitor signal before processing
    plot(TimePoints,SPL_mat,'k-');
    title ('Original input signal');
    ax=gca;
    xlabel('Time (s)');
    ylabel('SPL (dB)');
    ax.XLim = [0 TimePoints(n)];
    hold off;
    
    subplot(2,2,2)
    
    % Time-varying Fluctuation strength
    plot(TimePoints,FiT,'r-');
    title ('Time-varying Fluctuation');
    ax=gca;
    xlabel('Time (s)');
    ylabel('Fluctuation (vacil)');
    ax.XLim = [0 TimePoints(n)];
    hold off;
    
    subplot(2,2,4)
    
    % Time-averaged fluctuation as a function of ERB frequencies
    semilogx(ERBSfreq,Fi,'r-');
    title ('Specific Fluctuation');
    ax=gca;
    xlabel('Frequency (Hz)');
    ylabel('Fluctuation (vacil)');
    ax.XLim = [0 ERBSfreq(75)];
    hold off;
    
    subplot(2,2,3)
    % Time-varying Fluctuation strength spectrogram
    imagesc(TimePoints, 1:1:75, fi_mat);
    set(gca,'YDir','normal');
    ax=gca;
    axis tight;
    ax.Title.String = 'Fluctuation Strength (vacil)';
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'ERB_{channel}';
    hold off;
 
     % ******** AARAE RESULTS LEAVES **********
    
    % Time-varying SPL results leaf
    doresultleaf(SPL_mat','SPL [dB]',{'time'},...
        'Time',TimePoints','s',true,...
        'fluctuationstype', {'SPL over time'}, 'categorical',[],...
        'name','Time_varying_SPL');
    
    % Time-varying fluctuation results leaf
    doresultleaf(FiT','Fluctuation [vacil]',{'time'},...
        'Time',TimePoints','s',true,...
        'fluctuationstype', {'Fluctuation over time'}, 'categorical',[],...
        'name','Time_varying_Fluctuation');
    
    % Time-averaged Specific fluctuation results leaf
    doresultleaf(Fi,'Specific Fluctuation [vacil]',{'ERB'},...
        'Frequency',ERBSfreq,'Hz',true,...
        'fluctuationstype', {'Fluctuation over frequency'}, 'categorical', [],...
        'name','Time_averaged_specific_fluctuation');
    
    % Specific fluctuation spectrogram results leaf
        doresultleaf(fi_mat','Specific Fluctuation [vacil/ERB]',{'time','Channel'},...
            'Channel',[1:1:75]','ERB',true,...
            'Time',TimePoints','s',true,...
            'fluctuationstype', {'Fluctuation over ERB channel'}, 'categorical', [],...
            'name','Time_varying_specific_fluctuation');

    OUT.funcallback.name = 'Fluctuation.m';
    OUT.funcallback.inarg = {fs,cal};
else
    % output numbers only
    OUT = F; % total fluctuation
    varargout{1} = SPL_mat; % time-varying SPL of original input signal
    varargout{2} = FiT; % specific fluctuation over time
    varargout{3} = Fi; % specific fluctuation over ERB channel
end


else
    OUT = [];
end

%**************************************************************************
% Copyright (c) <2015>, <Ella Manor>
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the <ORGANISATION> nor the names of its contributors
%    may be used to endorse or promote products derived from this software 
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%**************************************************************************

