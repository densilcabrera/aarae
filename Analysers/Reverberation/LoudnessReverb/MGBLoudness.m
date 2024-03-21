function  [InstantaneousLoudness, ShortTermLoudness, LongTermLoudness, times] = MGBLoudness(signal,fs,filtermethod,level,Ltype, faster,decay,doplot)
% by Densil Cabrera & Doheon Lee (2009-10)
% Loudness calculated following Moore,Glasberg & Baer (1997)
% with Glasberg & Moore (2002) dynamic loudness
% and Moore & Glasberg (2007) binaural loudness
%
% WARNING 1 - this function is slow. Test it on a short duration sound file
% first, such as 0.5 seconds.
% However, the function is more efficiet than Densil's old version due to:
% * removing unnecessary calculations from the main program loop, and
% instead precalculating values before the loop. This uses much more memory;
% * vectorising all of the old for loops (except the time-step loop). This
% also uses much more memory.
% The function is also potentially faster than a straight implementation of
% Glasberg and Moore (2002) due to combining closely-spaced high frequency spectral components
% (see "faster" below). This is similar to the 'compact spectrum' used in
% Densil's old code (but implemented more efficiently, and with less data reduction).
% Another speed-up is to caclulate the FFT for low frequency components
% less often than for high frequency components. This, and the spectrum
% compression mentioned above, are both done when the 'faster' option is
% chosen.
%
% WARNING 2 - this code has tested well to the extent that we have
% validated it, but you use it at your own risk. We suggest that you test it for
% yourself.
%
% ISSUES
% We have included the 15.625 Hz component from the FFT in the loudness
% calculation, as it is affected by sound above 20 Hz. We found that our
% results matched the published results better when very low frequency
% tones were tested with this component included.
%
% We also obtained a closer match to published results, in the very low
% frequency range, when we used outer/middle ear transfer functions with
% fir2 filter values specified to match the published charts (rather than
% designing the measured transfer function to match the published charts).
% Both options are available through the selection of 'filtermethod'.
%
% The spectrum analysis procedure (of multiple FFTs) described by Glasberg
% and Moore (2002) causes some spectral spreading in the higher frequency range
% due to the zero-padded window functions. The program accounts for this,
% so that a 1 kHz 40 dB tone yields 1 sone, but this introduces a discrepancy
% of about 0.89 dB in the values that enter the loudness calculation.
% The discrepancy is probably negligible for most applications. An
% alternative approach to calibration is to add 0.89 dB relative to the rms
% value measured over the entire input wave.
%
% The time domain filters are probably preferable to the spectrum intensity weighting
% filters. The time domain filters introduce smoothing into
% transient parts of the wave (for example, a stopped sinusoid), which
% results in a smoother short-term loudness function.
%
% In this version, the modifications from Glasberg and Moore (2006) are not
% implemented, except for the revised middle ear transfer function (which
% is implemented as filtermethod 3).
%
% INPUTS
% "signal" is audio data (1 channel as a single column, or 2 channels as 2
% columns). If the data are 1 channel, then the calculation is for diotic
% loudness, as per Moore, Glasberg & Baer (1997). If the data are 2
% channels, then the calculation is for binaural loudness, as per Moore &
% Glasberg (2007).
%
% "fs" is the audio sampling rate. Note that this function resamples the
% audio to fs = 32 kHz, as this appears to be the sampling rate used by
% Glasberg and Moore (2002).
%
% "filtermethod" is the method of outer+middle ear, or only middle ear,
% filtering that is used. The methods in the current code are:
% 1 - combined outer and middle ear
% implemented as FIR2, based on the description and chart in Glasberg
% and Moore (2002). The measured transfer function of the filter is a good
% match to the chart.
% 2 - same as 1, except that the values input to fir2 match the chart,
% rather than the measured transfer function. This seems to produce a
% better match with the published results.
% 3 - middle ear only
% implemented as FIR2, based on the chart in Moore, Glasberg and Baer
% (1997).
% 4 - middle ear only, based on the chart in Glasberg and Moore (2006).
% 11 - combined outer and middle ear
% implemented by multiplying intensities in the frequency domain. Note that
% this has no temporal smearing because it is implemented after the fft.
% 12 - middle ear only, based on the 1997 chart,
% implemented by multiplying intensities in the frequency domain. Note that
% this has no temporal smearing because it is implemented after the fft.
% 13 - middle ear only, based on the 2006 chart.
% implemented by multiplying intensities in the frequency domain. Note that
% this has no temporal smearing because it is implemented after the fft.
%
% "level" is the SPL of a desired listening level.
%
% "Ltype" is a type of LEVEL. Options are "Leq" or "LAeq". 
% 
% "faster" uses 346 spectrum components (instead of 960), which speeds up the
% program considerably. This also reduces the memory requirement.
% The reduced number of components is implemented by maintaining a 0.1 Erb
% or less spacing of components in the high frequency range. Intensities
% are summed to form the new, less densely spaced, components. This may
% introduce some random error because the phase relationship between
% components does have the potential to influence the true sum (especially
% considering the zero-padding of the fft window at high frequencies).
% However, the Moore, Glasberg and Baer loudness model is not phase-
% sensitive, so this approximation is similar to that already inherent in
% the loudness model. Note that the 346-component spectrum is
% finer resolution than the 108-component 'compact spectrum' used in
% Densil's old implementation of Moore, Glasberg & Baer (1997).
% Use a value of 1 for "faster" if you wish to use it, otherwise use
% a value of 0 or omit it from the input arguments.
%
% "faster" also speeds up the calculation by conducting the FFT less
% frequently for low frequencies. The lowest part of the spectrum (15-80 Hz)
% is calculated once every 32 ms; the next part (80-500 Hz) every 16 ms, 
% and so on. The result is a somewhat stepped or
% jagged instantaneous loudness function when there are transients.
% However, it appears to have minimal effect on short term loudness.
%
% Another approach to speed-up might be to discard low intensity spectral
% components. However, to do this without using logical indexing would
% involve a significant redesign of the program structure, and would
% sacrifice the greater speed obtained by pre-calculating values that are
% not signal-dependant. As logical indexing introduces some delay,
% discarding low intensity values through logical indexing is
% unlikely to result in a net speed-up, except perhaps for analysing narrow
% band signals such as pure tones. The approach taken in this
% implementation is well-suited to broad-band signal analysis.
%
% "decay" is the time in milliseconds of additional decay after the end of
% the audio wave. During this period it is assumed that instantaneous
% loudness is 0, but short term and long term loudness are calculated based
% on their exponential decay constants. Note that 63 ms of silence is
% always added to the end of the file for the full loudness calculation,
% and "decay" adds silence in addition to this.
%
% "doplot" specifies whether a plot is produced at the end.
% 0 (or ommiting doplot) gives no plot
% 1 gives a plot of loudness (in sones)
% 2 gives a plot of loudness level (in phons)
% 3 does the same plot as 2, and also returns the loudness level instead of
% loudness in the function's output arguments
%
% OUTPUTS
% Output is instantaneous loudness, short term loudness and long term
% loudness, as described in Glasberg and Moore (2002). Times can also be output,
% although these are simple to calculate without using the 'times' output. Note that in this
% implementation, windows are centred on the indicated time
% (given by 'times', the units for which are seconds). The first
% window (0 ms) includes half a window of zero-padding before the start of
% the input wave.
%
% EXAMPLE OF CALLING THIS FUNCTION
% if you import a .wav file, you will have 'data' and 'fs'.
% Then you can analyse using:
% [InstantaneousLoudness, ShortTermLoudness, LongTermLoudness,times] =MGBLoudness(data,fs,1,0,1,1000,1);
% The example above uses time domain filtering for outer and middle
% ear transfer function, 0 dB calibration adjustment, the "faster"
% spectrum method, adds 1000 ms of decay time after the wave file,
% and produces a loudness plot.
%
% KEY SOURCES:
%
% B.R. Glasberg and B.C.J. Moore. 1990.
% Derivation of Auditory Filter Shapes from Notched Noise Data
% Hearing Research, 47: 103-137.
%
% B.C.J. Moore, B.R. Glasberg and T. Baer. 1997.
% A Model for the Prediction of Thresholds, Loudness, and Partial Loudness
% Journal of the Audio Engineering Society, 45(4): 224-240.
%
% B.R. Glasberg and B.C.J. Moore. 2002.
% A Model of Loudness Applicable to Time-Varying Sounds
% Journal of the Audio Engineering Society, 50(5): 331-342.
%
% B.R. Glasberg and B.C.J. Moore. 2006.
% Prediction of absolute thresholds and equal-loudness contours using a
% modified loudness model
% J. Acoust. Soc. Am. 120: 585-588
%
% B.C.J. Moore and B.R. Glasberg. 2007
% Modelling Binaural Loudness
% J. Acoust. Soc. Am. 121: 1604-1612

if nargin<8, doplot=0; end
if nargin<7, decay=500; end
if nargin<6, faster=1; end
if nargin<5, Ltype='LAFmax'; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(signal)
    [signal, fs]=audioread(char(signal));
end

% resample to fs = 32 kHz because Glasberg & Moore appear to be using this
% sampling rate in their paper. Please note that some of the calculations
% later in this function are hard-coded for fs = 32 kHz.
if fs ~= 32000
    signal = resample(signal,32000,fs);
    fs = 32000;
end

% check data
[len, nchan] = size(signal);
if nchan > len
    signal = signal';
    [len, nchan] = size(signal);
end
if nchan > 1
    signal = signal(:,1:2);
    nchan = 2;
end

% Level adjustment
if strcmpi(Ltype, 'Leq')==1;
    cal=level-20*log10(rms(signal(:)));
    
elseif strcmpi(Ltype, 'LAeq')==1;

    % A-weighting filter design
    % h=fdesign.audioweighting('WT,Class', 'A', 1, fs);
    % hd=design(h);
    weightType = 'A-weighting';
    weightFilt = weightingFilter(weightType,fs);

    % Fast integration (125 ms) filter design
    E = exp(-1/(0.125*fs)); % exponential term
    bb = 1 - E; % filter numerator (adjusts gain to compensate for denominator)
    a = [1 -E];

    % Applying the A-weighting filter & fast integration filter
    %LA=filter(hd, signal); % apply A-weigting
    LA = weightFilt(signal);
    LAF=filter(bb,a,abs(LA)); % rectify and apply temporal integration


    
    % Calculating the caliberation level (Cal)
    LAFmax=20*log10(mean(max(abs(LAF)))); % (mean maxima for 2-chan)
    cal=level-LAFmax;
else
    error('Ltype must be either Leq or LAeq');
end

% signal calibration offset
cal = cal-59.9259;
if length(level) == 1
    signal = signal .* 10.^(cal/20);
else
    signal(:,1) = signal(:,1) .* 10.^(cal(1)/20);
    signal(:,2) = signal(:,2) .* 10.^(cal(2)/20);
end

% use if you wish to monitor the sound pressure level here
% disp(['rms level of the entire wave ', num2str(10*log10(mean(signal.^2)+10e-99)+59.9259), ' dB'])
% the offset of 59.9259 dB was chosen so that a 1 kHz tone at 40 dB yields
% 1 sone. However, arguably 0.89 dB should be added to this (see comments
% elsewhere).   

% zeropad signal by 32 ms at start and 63 ms at end
signal = [zeros(1024,nchan);signal;zeros(2016,nchan)];

% outer and middle ear transfer function
% implemented as FIR filters
% alternatives to these can be used with frequency domain intensity
% scaling
switch filtermethod
    case 1
        % Free field 0 deg, outer and middle ear, from Fig 1 of Glasberg &
        % Moore 2002.
        % This is the filter method described in Glasberg and Moore 2002
        % Note that the values below do not match the chart in the very low
        % frequency range, but the measured response of the filter does
        % match the chart.
        % Filter refined by Doheon Lee
        %
        filtHzdB = [0 -200;
            20.00 -45.4;
            35.35 -28.06;
            49.60 -19.25;
            98.27 -12.68;
            393.81 -3.09;
            620.72 -0.44;
            744.77 5.14e-2;
            1.03e3 9.3e-2;
            1.24e3 0.58;
            1.57e3 2.02;
            1.99e3 3.6;
            3.11e3 8.58;
            3.99e3 8.63;
            4.99e3 5.13;
            5.95e3 -0.3;
            6.95e3 -4.9;
            7.94e3 -9.31;
            8.98e3 -11.38;
            9.96e3 -11.38;
            1.13e4 -9.83;
            1.24e4 -8.96;
            1.4e4 -14.58;
            1.49e4 -19.09;
            1.56e4 -16.22];

        b = fir2(4097,[filtHzdB(:,1)./(0.5*fs);1],[10.^(filtHzdB(:,2)./20);0]);
        signal = filter(b,1,signal);
        
    case 2
        % same as filtermethod 1, except that the specified values match the chart
        % in the low frequency range. However the measured response does not yield the profile
        % of the filter response with the transfer function (middle ear
        % transfer function combined with the outer ear transfer function)
        % in Glasber and Moore 2002
        
        filtHzdB = [0 -200;
            20.00 -39.4;
            35.35 -25.06;
            49.60 -19.25;
            98.27 -12.68;
            393.81 -3.09;
            620.72 -0.44;
            744.77 5.14e-2;
            1.03e3 9.3e-2;
            1.24e3 0.58;
            1.57e3 2.02;
            1.99e3 3.6;
            3.11e3 8.58;
            3.99e3 8.63;
            4.99e3 5.13;
            5.95e3 -0.3;
            6.95e3 -4.9;
            7.94e3 -9.31;
            8.98e3 -11.38;
            9.96e3 -11.38;
            1.13e4 -9.83;
            1.24e4 -8.96;
            1.4e4 -14.58;
            1.49e4 -19.09;
            1.56e4 -16.22];

        b = fir2(4097,[filtHzdB(:,1)./(0.5*fs);1],[10.^(filtHzdB(:,2)./20);0]);
        signal = filter(b,1,signal);
        
    case 3
        % just use middle ear transfer function (assume that the outer ear
        % transfer function has already been applied to the signal).
        % The values below were obtained by visual reading of Figure 3
        % in the paper by Moore, Glasberg and Baer (1997).
        % Since the figure in the paper represents the middle ear effective
        % attenuation, the data below is obtained by the effective attenuation
        % multiplied by -1.
        % Values were derived by Doheon Lee
        %
        filtHzdB = [0 -200;
            20.00 -50.2;
            22.20 -46.6;
            22.21 -46.16;
            23.38 -44.38;
            25.25 -41.04;
            25.5 -34.38;
            33.67 -23.8;
            36.56 -22.13;
            37.32 -21.54;
            41.38 -20.12;
            47.81 -18.31;
            51.4 -17.51;
            84.45 -13.77;
            101.73 -12.17;
            258.19 -6.9;
            339.7 -5.63;
            415.64 -4.35;
            485.48 -3.71;
            645.48 -2.91
            823.48 -2.58;
            1.25e3 -2.77;
            1.38e3 -3.25;
            1.51e3 -3.93;
            2.01e3 -8.72;
            2.48e3 -10.71;
            2.72e3 -9.62;
            2.92e3 -7.48;
            3.11e3 -6.71;
            3.34e3 -6.03;
            3.69e3 -5.7;
            4.8e3 -5.71;
            5.19e3 -5.95;
            5.67e3 -6.9;
            8.04e3 -11.27;
            1.0e4 -9.97;
            1.51e4 -17.31;
            1.60e4 -17.67];

        b = fir2(4096,filtHzdB(:,1)./(0.5*fs),10.^(filtHzdB(:,2)./20));
        signal = filter(b,1,signal);

    case 4 % Modified middle ear function (Glasberg and Moore 2006) - refined by Doheon Lee
        filtHzdB = [0 -120;
            20.0 -70.41;
            30.43 -28.6;
            36.56 -24.0;
            41.63 -21.77;
            45.98 -19.46;
            54.83 -17.49;
            73.35 -14.67;
            158.96 -9.6;
            217.6 -7.91;
            304.81 -6.33;
            417.28 -4.56;
            486.35 -3.77;
            707.98 -2.99;
            884.16 -2.46;
            1.02e3 -2.66;
            1.18e3 -3.97;
            1.36e3 -5.09;
            1.52e3 -5.81;
            2.47e3 -10.27;
            2.86e3 -8.83;
            2.99e3 -7.19;
            3.65e3 -6.53;
            4.46e3 -6.73;
            5.2e3 -7.51;
            6.06e3 -9.55;
            7.01e3 -11.78;
            7.81e3 -12.44;
            9.0e3 -10.66;
            9.82e3 -10.07;
            1.08e4 -12.04;
            1.17e4 -13.85;
            1.21e4 -14.63;
            1.25e4 -15.13;
            1.60e4 -20];

        b = fir2(4096,filtHzdB(:,1)./(0.5*fs),10.^(filtHzdB(:,2)./20));
        signal = filter(b,1,signal);

        % case 4
        % ADD MORE FILTERS HERE AS DESIRED
end % switch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PRECALCULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% precalculate values that are used in the main loop and which do not
% depend upon the input signal

% generate window functions for the six FFT bands
w1 = repmat(hann(2048),1,nchan);
w2 = repmat([zeros(512,1);hann(1024);zeros(512,1)],1,nchan);
w3 = repmat([zeros(768,1);hann(512);zeros(768,1)],1,nchan);
w4 = repmat([zeros(896,1);hann(256);zeros(896,1)],1,nchan);
w5 = repmat([zeros(960,1);hann(128);zeros(960,1)],1,nchan);
w6 = repmat([zeros(992,1);hann(64);zeros(992,1)],1,nchan);

f1 = ((2:961)-1)*fs/2048; % frequencies from FFTs in Hz
% nfreq1 = length(f1); % number of frequencies
% below are frequencies used in loudness calculation
if faster==0
    f = f1; % frequencies in Hz
else
    f = [15.625:15.625:2671.875,2695.3125:31.25:4132.8125,4171.875:46.875:5578.125,5632.8125:62.5:7007.8125,7078.125:78.125:8406.25,8492.1875:93.75:9898.4375,10000:109.375:11312.5,11429.6875:125:12804.6875,12937.5:140.625:14625,14773.4375,14929.6875]; % frequencies in Hz
    % 346 components < 0.1 erb spacing where possible
end
nfreq = length(f); % number of frequencies

% pre-allocate intensity array
I = zeros(length(f),nchan);
if faster ==1, IX = I; end
%I=zeros(nfreq1,1); 

% roex component levels
erb = 24.673 * (4.368 .* f./1000 + 1);
p = 4 .* f ./ erb;
g = abs((repmat(f',1,nfreq)-repmat(f,nfreq,1))./repmat(f,nfreq,1)); %memory hungry but fast
[row col] = find(g <= 2);
indg = find(g <= 2);
% row is "component" in the old code
% col is "i" in the old code

% excitation
p511k = 30.2012922; %p51 @ 1kHz - see Glasberg and Moore (1990)
ERBS = 2:0.25:39;
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
[row_g2_neg col_g2_neg]=find(g2<0);

% Binaural
B=0.08;
p_bi=1.5978;
erbx=repmat(ERBS, nERBS,1);
erby=repmat(ERBS',1, nERBS);
g_bi=(erbx-abs(erbx-erby))./erbx;
w_bi=exp(-B.*g_bi).^2;

% loudness
C = 0.047;
if nchan == 2
    C = C/0.75; % modified value from Moore and Glasberg 2007
end % if nchan == 2
Ethrq = ones(nERBS,1).*2.31;
alpha = ones(nERBS,1).*0.2;
A = ones(nERBS,1).*4.62;
g_N = ones(nERBS,1); % g
%note - 1:35 are characteristic frequencies less than 500 Hz
Ethrq(1:35) = 10.^(0.1.*((-4.500239673E-03.*ERBSfreq(1:35) + 3.666468615 + 1272.362339./ERBSfreq(1:35))));
g_N(1:35)=2.31./Ethrq(1:35);
g_NdB = 10*log10(g_N);
A(1:35) = -1.03703703E-04 .* (g_NdB(1:35).^3) -6.03174603E-04 .* (g_NdB(1:35).^2) ...
    -1.21375661E-01 .* g_NdB(1:35) + 4.58825396;
alpha(1:35) = 1.346860977E-23 .* (g_NdB(1:35).^3) +  2.571428571E-05 .* (g_NdB(1:35).^2) ...
    -2.071428571E-03 .* g_NdB(1:35) + 1.997142857E-01;
SpecLoudL = zeros(nERBS,1); SpecLoudR=zeros(nERBS,1);


% short-term loudness
XSa = 0.045; % attack constant
XSr = 0.02; % release constant
% long-term loudness
XLa = 0.01; % attack constant
XLr = 0.0005; % release constant

nwindows = floor(length(signal)./32 - 64); % number of windows
times = ((1:nwindows+2+decay)-1)'./1000; % time in seconds
% pre-allocate output
minloud = 0.003; % minimum loudness permitted
InstantaneousLoudness = ones(nwindows+2,1) .* minloud;
ShortTermLoudness = ones(nwindows+2,1) .* minloud;
LongTermLoudness = ones(nwindows+2,1) .* minloud;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spectrum analysis as described in Glasberg and Moore (2002)
for windowcount = 1:nwindows

    % acquire window of data
    x = signal((windowcount-1)*32+1:(windowcount-1)*32+2048,:);
    
    
    % Conduct FFTs and convert to intensities
    
    if faster==1
    %4050 - 15000 Hz
    X6 = fft(x .* w6);
    IX(260:960,:) = 32 * abs(X6(261:961,:)).^2;
    
        if rem(windowcount, 2) == 1
            %2540 - 4050 Hz
            X5 = fft(x .* w5);
            IX(163:259,:) = 16 * abs(X5(164:260,:)).^2;
            I(163:171,:) = IX(163:171,:); % 'low' freq part of compressed spectrum - see below
            
            if rem(windowcount, 4) == 1
                % 1250 - 2540 Hz
                X4 = fft(x .* w4);
                I(80:162,:) = 8 * abs(X4(81:163,:)).^2;

                if rem(windowcount, 8) == 1
                    % 500 - 1250 Hz
                    X3 = fft(x .* w3);
                    I(32:79,:) = 4 * abs(X3(33:80,:)).^2;

                    if rem(windowcount, 16) == 1
                        % 80 - 500 Hz
                        X2 = fft(x .* w2);
                        I(6:31,:) = 2 * abs(X2(7:32,:)).^2;

                        if rem(windowcount, 32) == 1
                            %20 - 80 Hz
                            X1 = fft(x .* w1);
                            I(1:5,:) = abs(X1(2:6,:)).^2;

                        end % if rem(windowcount, 32)
                    end % if rem(windowcount, 16)
                end % if rem(windowcount, 8)
            end % if rem(windowcount, 4)
        end % if rem(windowcount, 2)

        % compress spectrum to 345 elements by summing intensities,
        % maintain <0.1 erb spacing where possible
        % I(1:170,:) = IX(1:170,:); this line is above for speed
        I(171:217,:)= sum(cat(3,IX(172:2:264,:), IX(173:2:265,:)),3);
        I(218:248,:)=sum(cat(3,IX(266:3:356,:), IX(267:3:357,:), IX(268:3:358,:)),3);
        I(249:271,:)=sum(cat(3,IX(359:4:447,:), IX(360:4:448,:), IX(361:4:449,:), IX(362:4:450,:)),3);
        I(272:289,:)=sum(cat(3,IX(451:5:536,:), IX(452:5:537,:), IX(453:5:538,:), IX(454:5:539,:), IX(455:5:540,:)),3);
        I(290:305,:)=sum(cat(3,IX(541:6:631,:), IX(542:6:632,:), IX(543:6:633,:), IX(544:6:634,:), IX(545:6:635,:), IX(546:6:636,:)),3);
        I(306:318,:)=sum(cat(3,IX(637:7:721,:), IX(638:7:722,:), IX(639:7:723,:), IX(640:7:724,:), IX(641:7:725,:), IX(642:7:726,:), IX(643:7:727,:)),3);
        I(319:330,:)=sum(cat(3,IX(728:8:816,:), IX(729:8:817,:), IX(730:8:818,:), IX(731:8:819,:), IX(732:8:820,:), IX(733:8:821,:), IX(734:8:822,:), IX(735:8:823,:)),3);
        I(331:343,:)=sum(cat(3,IX(824:9:932,:), IX(825:9:933,:), IX(826:9:934,:), IX(827:9:935,:), IX(828:9:936,:), IX(829:9:937,:), IX(830:9:938,:), IX(831:9:939,:), IX(832:9:940,:)),3);
        I(344,:)=sum(IX(941:950,:),1);
        I(345,:)=sum(IX(951:960,:),1);
        
    else % if ~faster
        %4050 - 15000 Hz
        X6 = fft(x .* w6);
        I(259:959,:) = 32 * abs(X6(261:961,:)).^2;
        %2540 - 4050 Hz
        X5 = fft(x .* w5);
        I(162:258,:) = 16 * abs(X5(164:260,:)).^2;
        % 1250 - 2540 Hz
        X4 = fft(x .* w4);
        I(79:161,:) = 8 * abs(X4(81:163,:)).^2;
        % 500 - 1250 Hz
        X3 = fft(x .* w3);
        I(31:78,:) = 4 * abs(X3(33:80,:)).^2;
        % 80 - 500 Hz
        X2 = fft(x .* w2);
        I(5:30,:) = 2 * abs(X2(7:32,:)).^2;
        %20 - 80 Hz
        X1 = fft(x .* w1);
        I(1:4,:) = abs(X1(3:6,:)).^2;
    end % if faster

    % 40 dB 1 kHz pure tone should
    % yield 1 sone -  for testing purposes
    % call the function using "faster" , 1 chan, if you wish to test this
    % the line below overwrites the "faster" intensity spectrum with a 1 kHz
    % tone at 40 dB
    %   I=[zeros(63,1);10^(40/10);zeros(282,1)];
    
%     To test 31.25 Hz at 100 dB (without Middle and outer ear trasnfer
%     function), you can use the following line

%     I=[zeros(1,1);10^(100/10);zeros(282+62,1)];

    % use if you wish to monitor the sound pressure level here
    %     disp(['level after windowing ', num2str(10*log10(sum(I)+1e-99)), ' dB'])
    % note - A 1 kHz pure tone at 39.13 dB yields a loudness of 1 sone
    % and this 0.89 dB discrepancy is due to the spectrum analysis
    % specified by Glasberg and Moore (2002). One approach to this would be
    % to use the following level instead:
    % disp(['level after windowing ', num2str(10*log10(sum(I)+1e-99)+0.89), ' dB'])
    
    %hold on
    %plot(f,10*log10(I(:,1)))
    

    switch filtermethod
        % The following filter methods are frequency-domain weighting of the intensity spectrum
        % this might be faster than the time domain method (especially for the compressed spectrum),
        % but will result in no temporal smearing as it acts on a single
        % spectrum window.
        % Note that Glasberg and Moore (2002) use time-domain filtering
        % (filtermethod=1), not the methods below.
        % Values are hard-coded for speed.
        % These filters are not significantly faster than the time domain
        % alternatives
        case 11
            % free field approximation combined with middle ear, from
            % chart in Glasberg & Moore (2002)
            if faster==0
                % Values refined by Doheon Lee
                I=I.*[0.0012912,0.0092023,0.017748,0.028844,0.04688,0.058618,0.065877,0.074035,0.083203,0.093506,0.10509,0.1181,0.13272,0.14916,0.16763,0.18839,0.21172,0.23793,0.2674,0.30051,0.33772,0.37954,0.42654,0.47936,0.50761,0.52939,0.55211,0.5758,0.60051,0.62628,0.65315,0.68118,0.71041,0.7409,0.77269,0.80585,0.84043,0.8765,0.90718,0.92021,0.93341,0.94681,0.9604,0.97419,0.98817,1.0024,1.0121,1.0126,1.0131,1.0137,1.0142,1.0147,1.0153,1.0158,1.0163,1.0169,1.0174,1.0179,1.0185,1.019,1.0195,1.0201,1.0206,1.0212,1.0223,1.0309,1.0395,1.0482,1.057,1.0659,1.0748,1.0838,1.0929,1.1021,1.1113,1.1206,1.13,1.1395,1.1544,1.1727,1.1912,1.2101,1.2292,1.2487,1.2684,1.2885,1.3089,1.3296,1.3507,1.372,1.3937,1.4158,1.4382,1.461,1.4841,1.5076,1.5314,1.5556,1.5803,1.6035,1.6253,1.6475,1.6699,1.6927,1.7157,1.7391,1.7628,1.7868,1.8112,1.8358,1.8609,1.8862,1.9119,1.938,1.9644,1.9912,2.0183,2.0458,2.0737,2.1019,2.1306,2.1596,2.189,2.2189,2.2491,2.2797,2.3144,2.3518,2.3897,2.4282,2.4674,2.5072,2.5476,2.5887,2.6304,2.6728,2.716,2.7597,2.8043,2.8495,2.8954,2.9421,2.9896,3.0378,3.0868,3.1365,3.1871,3.2385,3.2907,3.3438,3.3977,3.4525,3.5082,3.5648,3.6222,3.6807,3.74,3.8003,3.8616,3.9239,3.9872,4.0515,4.1168,4.1832,4.2506,4.3192,4.3888,4.4596,4.5315,4.6046,4.6788,4.7543,4.831,4.9089,4.988,5.0685,5.1502,5.2332,5.3176,5.4034,5.4905,5.5791,5.669,5.7605,5.8533,5.9477,6.0436,6.1411,6.2401,6.3408,6.443,6.5469,6.6525,6.7598,6.8688,6.9795,7.0921,7.2065,7.2125,7.214,7.2154,7.2169,7.2184,7.2199,7.2213,7.2228,7.2243,7.2258,7.2272,7.2287,7.2302,7.2317,7.2332,7.2346,7.2361,7.2376,7.2391,7.2406,7.242,7.2435,7.245,7.2465,7.248,7.2494,7.2509,7.2524,7.2539,7.2554,7.2569,7.2583,7.2598,7.2613,7.2628,7.2643,7.2658,7.2672,7.2687,7.2702,7.2717,7.2732,7.2747,7.2762,7.2777,7.2791,7.2806,7.2821,7.2836,7.2851,7.2866,7.2881,7.2896,7.2911,7.2925,7.294,7.236,7.1455,7.0561,6.9678,6.8806,6.7945,6.7095,6.6255,6.5426,6.4607,6.3799,6.3,6.2212,6.1434,6.0665,5.9906,5.9156,5.8416,5.7685,5.6963,5.625,5.5546,5.4851,5.4165,5.3487,5.2818,5.2157,5.1504,5.086,5.0223,4.9595,4.8974,4.8362,4.7756,4.7159,4.6569,4.5986,4.541,4.4842,4.4281,4.3727,4.318,4.264,4.2106,4.1579,4.1059,4.0545,4.0038,3.9537,3.9042,3.8553,3.8071,3.7595,3.7124,3.666,3.6201,3.5748,3.5301,3.4859,3.4423,3.3992,3.3567,3.3146,3.2732,3.2162,3.1514,3.0879,3.0257,2.9648,2.9051,2.8465,2.7892,2.733,2.678,2.624,2.5711,2.5194,2.4686,2.4189,2.3701,2.3224,2.2756,2.2298,2.1849,2.1408,2.0977,2.0555,2.0141,1.9735,1.9337,1.8948,1.8566,1.8192,1.7826,1.7467,1.7115,1.677,1.6432,1.6101,1.5777,1.5459,1.5148,1.4842,1.4543,1.425,1.3963,1.3682,1.3406,1.3136,1.2872,1.2612,1.2358,1.2109,1.1866,1.1627,1.1392,1.1163,1.0938,1.0718,1.0502,1.029,1.0083,0.98798,0.96807,0.94857,0.93017,0.9149,0.89989,0.88512,0.87059,0.8563,0.84224,0.82842,0.81482,0.80145,0.78829,0.77535,0.76263,0.75011,0.7378,0.72569,0.71378,0.70206,0.69054,0.6792,0.66806,0.65709,0.6463,0.6357,0.62526,0.615,0.60491,0.59498,0.58521,0.57561,0.56616,0.55687,0.54772,0.53873,0.52989,0.52119,0.51264,0.50423,0.49595,0.48781,0.4798,0.47193,0.46418,0.45656,0.44907,0.4417,0.43445,0.42732,0.4203,0.4134,0.40662,0.39994,0.39338,0.38692,0.38057,0.37433,0.36818,0.36214,0.35619,0.35035,0.3446,0.33894,0.33338,0.32791,0.32256,0.31743,0.31238,0.30742,0.30253,0.29772,0.29299,0.28833,0.28374,0.27923,0.27479,0.27042,0.26612,0.26189,0.25773,0.25363,0.2496,0.24563,0.24173,0.23788,0.2341,0.23038,0.22672,0.22311,0.21956,0.21607,0.21264,0.20926,0.20593,0.20266,0.19943,0.19626,0.19314,0.19007,0.18705,0.18408,0.18115,0.17827,0.17544,0.17265,0.1699,0.1672,0.16454,0.16193,0.15935,0.15682,0.15433,0.15187,0.14946,0.14708,0.14474,0.14244,0.14018,0.13795,0.13575,0.1336,0.13147,0.12938,0.12733,0.1253,0.12331,0.12135,0.11942,0.11752,0.11652,0.11569,0.11486,0.11404,0.11323,0.11242,0.11162,0.11082,0.11003,0.10924,0.10846,0.10769,0.10692,0.10616,0.1054,0.10465,0.1039,0.10316,0.10243,0.10169,0.10097,0.10025,0.099533,0.098823,0.098118,0.097418,0.096723,0.096033,0.095347,0.094667,0.093991,0.093321,0.092655,0.091994,0.091337,0.090686,0.090039,0.089396,0.088758,0.088125,0.087496,0.086872,0.086252,0.085636,0.085025,0.084419,0.083816,0.083218,0.082624,0.082035,0.08145,0.080868,0.080291,0.079718,0.07915,0.078585,0.078024,0.077467,0.076915,0.076366,0.075821,0.07528,0.074743,0.074209,0.07368,0.073154,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072948,0.073252,0.073557,0.073864,0.074172,0.074482,0.074792,0.075104,0.075417,0.075732,0.076048,0.076365,0.076683,0.077003,0.077324,0.077647,0.07797,0.078296,0.078622,0.07895,0.079279,0.07961,0.079942,0.080275,0.08061,0.080946,0.081284,0.081623,0.081963,0.082305,0.082648,0.082993,0.083339,0.083686,0.084035,0.084386,0.084738,0.085091,0.085446,0.085802,0.08616,0.08652,0.08688,0.087243,0.087607,0.087972,0.088339,0.088707,0.089077,0.089449,0.089822,0.090196,0.090572,0.09095,0.091329,0.09171,0.092093,0.092477,0.092862,0.09325,0.093638,0.094029,0.094421,0.094815,0.09521,0.095607,0.096006,0.096406,0.096808,0.097212,0.097618,0.098025,0.098433,0.098844,0.099256,0.09967,0.10009,0.1005,0.10092,0.10134,0.10177,0.10219,0.10262,0.10304,0.10347,0.10391,0.10423,0.10453,0.10482,0.10512,0.10542,0.10572,0.10602,0.10633,0.10663,0.10693,0.10724,0.10754,0.10785,0.10816,0.10847,0.10877,0.10908,0.10939,0.10971,0.11002,0.11033,0.11065,0.11096,0.11128,0.1116,0.11191,0.11223,0.11255,0.11287,0.11319,0.11352,0.11384,0.11417,0.11449,0.11482,0.11514,0.11547,0.1158,0.11613,0.11646,0.11679,0.11713,0.11746,0.1178,0.11813,0.11847,0.11881,0.11914,0.11948,0.11982,0.12017,0.12051,0.12085,0.1212,0.12154,0.12189,0.12223,0.12258,0.12293,0.12328,0.12363,0.12399,0.12434,0.12469,0.12505,0.12541,0.12576,0.12612,0.12648,0.12684,0.12642,0.12483,0.12326,0.12171,0.12019,0.11868,0.11719,0.11571,0.11426,0.11283,0.11141,0.11001,0.10863,0.10726,0.10592,0.10459,0.10327,0.10198,0.1007,0.099432,0.098184,0.096951,0.095733,0.094531,0.093344,0.092172,0.091014,0.089871,0.088743,0.087628,0.086528,0.085441,0.084368,0.083309,0.082263,0.08123,0.08021,0.079202,0.078208,0.077226,0.076256,0.075298,0.074353,0.073419,0.072497,0.071587,0.070688,0.0698,0.068923,0.068058,0.067203,0.066359,0.065526,0.064703,0.063891,0.063088,0.062296,0.061514,0.060741,0.059979,0.059225,0.058482,0.057747,0.057022,0.056306,0.055599,0.054901,0.054211,0.053531,0.052858,0.052195,0.051539,0.050892,0.050253,0.049622,0.048999,0.048383,0.047776,0.047176,0.046583,0.045998,0.045421,0.04485,0.044287,0.043731,0.043182,0.04264,0.042104,0.041575,0.041053,0.040538,0.040029,0.039526,0.03903,0.03854,0.038056,0.037578,0.037106,0.03664,0.03618,0.035725,0.035277,0.034834,0.034211,0.0336,0.033,0.03241,0.031831,0.031262,0.030704,0.030155,0.029616,0.029087,0.028567,0.028057,0.027556,0.027063,0.02658,0.026105,0.025638,0.02518,0.02473,0.024289,0.023855,0.023428,0.02301,0.022599,0.022195,0.021798,0.021409,0.021026,0.020651,0.020282,0.019919,0.019563,0.019214,0.018871,0.018533,0.018202,0.017877,0.017558,0.017244,0.016936,0.016633,0.016336,0.016044,0.015757,0.015476,0.015199,0.014928,0.014661,0.014399,0.014142,0.013889,0.013641,0.013397,0.013158,0.012923,0.012692,0.012465,0.012404,0.012588,0.012775,0.012965,0.013158,0.013353,0.013552]';
            elseif faster==1
                I=I.*[0.0012912,0.0092023,0.017748,0.028844,0.04688,0.058618,0.065877,0.074035,0.083203,0.093506,0.10509,0.1181,0.13272,0.14916,0.16763,0.18839,0.21172,0.23793,0.2674,0.30051,0.33772,0.37954,0.42654,0.47936,0.50761,0.52939,0.55211,0.5758,0.60051,0.62628,0.65315,0.68118,0.71041,0.7409,0.77269,0.80585,0.84043,0.8765,0.90718,0.92021,0.93341,0.94681,0.9604,0.97419,0.98817,1.0024,1.0121,1.0126,1.0131,1.0137,1.0142,1.0147,1.0153,1.0158,1.0163,1.0169,1.0174,1.0179,1.0185,1.019,1.0195,1.0201,1.0206,1.0212,1.0223,1.0309,1.0395,1.0482,1.057,1.0659,1.0748,1.0838,1.0929,1.1021,1.1113,1.1206,1.13,1.1395,1.1544,1.1727,1.1912,1.2101,1.2292,1.2487,1.2684,1.2885,1.3089,1.3296,1.3507,1.372,1.3937,1.4158,1.4382,1.461,1.4841,1.5076,1.5314,1.5556,1.5803,1.6035,1.6253,1.6475,1.6699,1.6927,1.7157,1.7391,1.7628,1.7868,1.8112,1.8358,1.8609,1.8862,1.9119,1.938,1.9644,1.9912,2.0183,2.0458,2.0737,2.1019,2.1306,2.1596,2.189,2.2189,2.2491,2.2797,2.3144,2.3518,2.3897,2.4282,2.4674,2.5072,2.5476,2.5887,2.6304,2.6728,2.716,2.7597,2.8043,2.8495,2.8954,2.9421,2.9896,3.0378,3.0868,3.1365,3.1871,3.2385,3.2907,3.3438,3.3977,3.4525,3.5082,3.5648,3.6222,3.6807,3.74,3.8003,3.8616,3.9239,3.9872,4.0515,4.1168,4.1832,4.2506,4.3192,4.3888,4.4596,4.5315,4.6046,4.7164,4.8698,5.0281,5.1916,5.3603,5.5346,5.7146,5.9004,6.0922,6.2903,6.4948,6.7059,6.9239,7.149,7.2132,7.2162,7.2191,7.2221,7.225,7.228,7.2309,7.2339,7.2369,7.2398,7.2428,7.2457,7.2487,7.2517,7.2546,7.2576,7.2606,7.2635,7.2665,7.2695,7.2725,7.2754,7.2784,7.2814,7.2844,7.2873,7.2903,7.2933,7.1906,7.0118,6.8374,6.6673,6.5015,6.3,6.0665,5.8416,5.625,5.4165,5.2157,5.0223,4.8362,4.6569,4.4842,4.318,4.1579,4.0038,3.8553,3.7124,3.5748,3.4423,3.3146,3.1514,2.9648,2.7892,2.624,2.4686,2.3224,2.1849,2.0555,1.9337,1.8192,1.7115,1.6101,1.5148,1.4106,1.3003,1.1987,1.105,1.0186,0.93897,0.87782,0.82159,0.76896,0.71971,0.67361,0.63046,0.59007,0.55228,0.5169,0.48379,0.4528,0.4238,0.39665,0.37124,0.34746,0.3252,0.30496,0.28374,0.26189,0.24173,0.22311,0.20593,0.19007,0.17544,0.16193,0.14946,0.13795,0.12733,0.11752,0.11323,0.10924,0.1054,0.10169,0.098118,0.094667,0.091011,0.087183,0.083517,0.080004,0.07664,0.073417,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.072778,0.073557,0.075732,0.07797,0.080275,0.082648,0.085091,0.087607,0.090196,0.092862,0.095607,0.098433,0.10134,0.10423,0.10648,0.10893,0.11144,0.114,0.11663,0.11931,0.12206,0.12487,0.12404,0.11212,0.10134,0.091591,0.082263,0.073419,0.065526,0.058482,0.052195,0.046583,0.041575,0.037106,0.03241,0.027556,0.023428,0.019919,0.016936,0.01427,0.012682]';
            end

        case 12
            %just the middle ear transfer function from Moore
            %Glasberg and Baer (1997) derived as above
            % Values refined by Doheon Lee
            if faster==0
                I=I.*[0.0026613,0.013889,0.023692,0.035598,0.051181,0.06438,0.072674,0.082037,0.092606,0.10454,0.118,0.13321,0.15037,0.16974,0.19161,0.20969,0.22178,0.23457,0.2481,0.2624,0.27786,0.29523,0.31369,0.3333,0.35414,0.37215,0.38462,0.39751,0.41084,0.42461,0.43277,0.44063,0.44863,0.45677,0.46506,0.4735,0.4821,0.49085,0.49976,0.50883,0.51404,0.51748,0.52094,0.52443,0.52794,0.53147,0.53503,0.53861,0.54221,0.54584,0.5495,0.55181,0.55093,0.55005,0.54917,0.54829,0.54741,0.54653,0.54566,0.54478,0.54391,0.54304,0.54217,0.5413,0.54044,0.53957,0.53871,0.53784,0.53698,0.53612,0.53526,0.53441,0.53355,0.5327,0.53184,0.53099,0.53014,0.52929,0.52845,0.52147,0.51459,0.5078,0.5011,0.49449,0.48796,0.48152,0.47517,0.46713,0.45843,0.44988,0.44149,0.43326,0.42518,0.41726,0.40948,0.39959,0.38605,0.37297,0.36033,0.34813,0.33633,0.32494,0.31393,0.30329,0.29302,0.28309,0.2735,0.26423,0.25528,0.24663,0.23828,0.2302,0.2224,0.21487,0.20759,0.20056,0.19376,0.1872,0.18085,0.17473,0.16881,0.16309,0.15756,0.15223,0.14707,0.14209,0.13727,0.13354,0.13152,0.12953,0.12758,0.12565,0.12375,0.12188,0.12004,0.11822,0.11643,0.11467,0.11294,0.11123,0.10955,0.10789,0.10626,0.10466,0.10307,0.10152,0.099982,0.09847,0.096981,0.095515,0.094071,0.092649,0.091249,0.089869,0.08851,0.087172,0.085855,0.085307,0.086713,0.088141,0.089593,0.091069,0.09257,0.094095,0.095645,0.09722,0.098822,0.10045,0.1021,0.10379,0.1055,0.10723,0.109,0.11308,0.11752,0.12213,0.12692,0.1319,0.13708,0.14246,0.14805,0.15386,0.1599,0.16618,0.1727,0.17896,0.18159,0.18426,0.18696,0.18971,0.1925,0.19532,0.19819,0.2011,0.20406,0.20705,0.21009,0.21318,0.21549,0.2178,0.22013,0.22248,0.22486,0.22727,0.2297,0.23215,0.23463,0.23714,0.23968,0.24224,0.24483,0.24745,0.24966,0.25051,0.25136,0.25222,0.25307,0.25393,0.2548,0.25566,0.25653,0.2574,0.25828,0.25915,0.26004,0.26092,0.26181,0.26269,0.26359,0.26448,0.26538,0.26628,0.26719,0.2681,0.26901,0.26915,0.26914,0.26913,0.26912,0.26911,0.2691,0.26909,0.26909,0.26908,0.26907,0.26906,0.26905,0.26904,0.26903,0.26902,0.26902,0.26901,0.269,0.26899,0.26898,0.26897,0.26896,0.26895,0.26895,0.26894,0.26893,0.26892,0.26891,0.2689,0.26889,0.26888,0.26888,0.26887,0.26886,0.26885,0.26884,0.26883,0.26882,0.26881,0.26881,0.2688,0.26879,0.26878,0.26877,0.26876,0.26875,0.26875,0.26874,0.26873,0.26872,0.26871,0.2687,0.26869,0.26868,0.26868,0.26867,0.26866,0.26865,0.26864,0.26863,0.26862,0.26861,0.26861,0.2686,0.26859,0.26858,0.26857,0.26856,0.26855,0.26854,0.26854,0.26806,0.26747,0.26687,0.26628,0.2657,0.26511,0.26452,0.26394,0.26335,0.26277,0.26219,0.26161,0.26103,0.26045,0.25988,0.2593,0.25873,0.25816,0.25759,0.25702,0.25645,0.25588,0.25532,0.25475,0.25419,0.25258,0.25079,0.24901,0.24724,0.24549,0.24375,0.24202,0.2403,0.2386,0.2369,0.23522,0.23355,0.2319,0.23025,0.22862,0.22699,0.22538,0.22378,0.2222,0.22062,0.21906,0.2175,0.21596,0.21443,0.2129,0.21139,0.20989,0.2084,0.20693,0.20546,0.20401,0.20266,0.20132,0.19999,0.19867,0.19736,0.19605,0.19475,0.19347,0.19219,0.19092,0.18965,0.1884,0.18715,0.18592,0.18469,0.18347,0.18225,0.18105,0.17985,0.17866,0.17748,0.17631,0.17514,0.17398,0.17283,0.17169,0.17056,0.16943,0.16831,0.16719,0.16609,0.16499,0.1639,0.16282,0.16174,0.16067,0.15961,0.15855,0.1575,0.15646,0.15543,0.1544,0.15338,0.15237,0.15136,0.15036,0.14936,0.14838,0.14739,0.14642,0.14545,0.14449,0.14353,0.14259,0.14164,0.14071,0.13978,0.13885,0.13793,0.13702,0.13612,0.13522,0.13432,0.13343,0.13255,0.13168,0.1308,0.12994,0.12908,0.12823,0.12738,0.12654,0.1257,0.12487,0.12404,0.12322,0.12241,0.1216,0.1208,0.12,0.1192,0.11842,0.11763,0.11685,0.11608,0.11531,0.11455,0.11379,0.11304,0.11229,0.11155,0.11081,0.11008,0.10935,0.10863,0.10791,0.1072,0.10649,0.10579,0.10509,0.10439,0.1037,0.10302,0.10233,0.10166,0.10099,0.10032,0.099655,0.098996,0.098342,0.097691,0.097045,0.096404,0.095766,0.095133,0.094504,0.093879,0.093259,0.092642,0.092029,0.091421,0.090816,0.090216,0.089619,0.089027,0.088438,0.087853,0.087273,0.086696,0.086122,0.085553,0.084987,0.084425,0.083867,0.083313,0.082762,0.082214,0.081671,0.081131,0.080594,0.080062,0.079532,0.079006,0.078484,0.077965,0.07745,0.076937,0.076429,0.075923,0.075421,0.074923,0.074723,0.074902,0.075081,0.07526,0.07544,0.07562,0.075801,0.075982,0.076163,0.076345,0.076528,0.076711,0.076894,0.077078,0.077262,0.077446,0.077631,0.077817,0.078003,0.078189,0.078376,0.078563,0.078751,0.078939,0.079128,0.079317,0.079506,0.079696,0.079887,0.080077,0.080269,0.080461,0.080653,0.080845,0.081039,0.081232,0.081426,0.081621,0.081816,0.082011,0.082207,0.082404,0.082601,0.082798,0.082996,0.083194,0.083393,0.083592,0.083792,0.083992,0.084193,0.084394,0.084595,0.084797,0.085,0.085203,0.085407,0.085611,0.085815,0.08602,0.086226,0.086432,0.086638,0.086845,0.087053,0.087261,0.087469,0.087678,0.087888,0.088098,0.088308,0.088519,0.088731,0.088943,0.089155,0.089368,0.089582,0.089796,0.09001,0.090225,0.090441,0.090657,0.090873,0.091091,0.091308,0.091526,0.091745,0.091964,0.092184,0.092404,0.092625,0.092846,0.093068,0.09329,0.093513,0.093737,0.093961,0.094185,0.09441,0.094636,0.094862,0.095088,0.095316,0.095543,0.095772,0.096,0.09623,0.09646,0.09669,0.096921,0.097153,0.097385,0.097617,0.097851,0.098084,0.098319,0.098554,0.098789,0.099025,0.099262,0.099499,0.099737,0.099975,0.10021,0.10045,0.10069,0.10017,0.099656,0.099141,0.098629,0.09812,0.097613,0.097109,0.096607,0.096108,0.095612,0.095118,0.094627,0.094138,0.093652,0.093168,0.092687,0.092208,0.091732,0.091258,0.090787,0.090318,0.089852,0.089388,0.088926,0.088467,0.08801,0.087555,0.087103,0.086653,0.086206,0.085761,0.085318,0.084877,0.084439,0.084003,0.083569,0.083137,0.082708,0.082281,0.081856,0.081433,0.081012,0.080594,0.080178,0.079764,0.079352,0.078942,0.078534,0.078129,0.077725,0.077324,0.076924,0.076527,0.076132,0.075739,0.075347,0.074958,0.074571,0.074186,0.073803,0.073422,0.073043,0.072665,0.07229,0.071917,0.071545,0.071176,0.070808,0.070442,0.070079,0.069717,0.069357,0.068998,0.068642,0.068288,0.067935,0.067584,0.067235,0.066888,0.066542,0.066199,0.065857,0.065517,0.065178,0.064842,0.064507,0.064174,0.063842,0.063512,0.063184,0.062858,0.062534,0.062211,0.061889,0.06157,0.061252,0.060935,0.060621,0.060307,0.059996,0.059686,0.059378,0.059071,0.058766,0.058463,0.058161,0.05786,0.057562,0.057264,0.056969,0.056674,0.056382,0.05609,0.055801,0.055513,0.055226,0.054941,0.054657,0.054375,0.054094,0.053814,0.053536,0.05326,0.052985,0.052711,0.052439,0.052168,0.051899,0.051631,0.051364,0.051099,0.050835,0.050572,0.050311,0.050051,0.049793,0.049536,0.04928,0.049025,0.048772,0.04852,0.04827,0.04802,0.047772,0.047526,0.04728,0.047036,0.046793,0.046551,0.046311,0.046072,0.045834,0.045597,0.045362,0.045127,0.044894,0.044662,0.044432,0.044202,0.043974,0.043747,0.043521,0.043296,0.043073,0.04285,0.042629,0.042409,0.04219,0.041972,0.041755,0.041539,0.041325,0.041111,0.040899,0.040688,0.040478,0.040269,0.040061,0.039854,0.039648,0.039443,0.039239,0.039037,0.038835,0.038635,0.038435,0.038237,0.038039,0.037843,0.037647,0.037453,0.037259,0.037067,0.036875,0.036685,0.036496,0.036307,0.03612,0.035933,0.035747,0.035563,0.035379,0.035196,0.035015,0.034834,0.034654,0.034475,0.034297,0.03412,0.033944,0.033768,0.033594,0.03342,0.033248,0.033076,0.032905,0.032735,0.032566,0.032398,0.032231,0.032064,0.031899,0.031734,0.03157,0.031407,0.031245,0.031083,0.030923,0.030763,0.030604,0.030446,0.030289,0.030133,0.029977,0.029822,0.029668,0.029515,0.029362,0.029211,0.02906,0.02891,0.028761,0.028612,0.028464,0.028317,0.028171,0.028025,0.027881,0.027737,0.027593,0.027451,0.027309,0.027168,0.027028,0.026888,0.026749,0.026611,0.026474,0.026337,0.026201,0.026066,0.025931,0.025797,0.025664,0.025531,0.0254,0.025268,0.025138,0.025008,0.024879,0.02475,0.024623,0.024495,0.024369,0.024243,0.024118,0.023993,0.023869,0.023746,0.023623,0.023501,0.02338,0.023259,0.023139,0.02302,0.022901,0.022783,0.022665,0.022548,0.022431,0.022315,0.0222,0.022086,0.021972,0.021858,0.021745,0.021633,0.021521,0.02141,0.021299,0.021189,0.02108,0.020971,0.020863,0.020755,0.020648,0.020541,0.020435,0.02033,0.020225,0.02012,0.020016,0.019913,0.01981,0.019708,0.019606,0.019505,0.019404,0.019304,0.019204]';
            elseif faster ==1
                % Below is the improved values
                I=I.*[0.0026613,0.013889,0.023692,0.035598,0.051181,0.06438,0.072674,0.082037,0.092606,0.10454,0.118,0.13321,0.15037,0.16974,0.19161,0.20969,0.22178,0.23457,0.2481,0.2624,0.27786,0.29523,0.31369,0.3333,0.35414,0.37215,0.38462,0.39751,0.41084,0.42461,0.43277,0.44063,0.44863,0.45677,0.46506,0.4735,0.4821,0.49085,0.49976,0.50883,0.51404,0.51748,0.52094,0.52443,0.52794,0.53147,0.53503,0.53861,0.54221,0.54584,0.5495,0.55181,0.55093,0.55005,0.54917,0.54829,0.54741,0.54653,0.54566,0.54478,0.54391,0.54304,0.54217,0.5413,0.54044,0.53957,0.53871,0.53784,0.53698,0.53612,0.53526,0.53441,0.53355,0.5327,0.53184,0.53099,0.53014,0.52929,0.52845,0.52147,0.51459,0.5078,0.5011,0.49449,0.48796,0.48152,0.47517,0.46713,0.45843,0.44988,0.44149,0.43326,0.42518,0.41726,0.40948,0.39959,0.38605,0.37297,0.36033,0.34813,0.33633,0.32494,0.31393,0.30329,0.29302,0.28309,0.2735,0.26423,0.25528,0.24663,0.23828,0.2302,0.2224,0.21487,0.20759,0.20056,0.19376,0.1872,0.18085,0.17473,0.16881,0.16309,0.15756,0.15223,0.14707,0.14209,0.13727,0.13354,0.13152,0.12953,0.12758,0.12565,0.12375,0.12188,0.12004,0.11822,0.11643,0.11467,0.11294,0.11123,0.10955,0.10789,0.10626,0.10466,0.10307,0.10152,0.099982,0.09847,0.096981,0.095515,0.094071,0.092649,0.091249,0.089869,0.08851,0.087172,0.085855,0.085307,0.086713,0.088141,0.089593,0.091069,0.09257,0.094095,0.095645,0.09722,0.098822,0.10045,0.1021,0.10379,0.10636,0.11092,0.1198,0.12939,0.13974,0.15093,0.16301,0.17605,0.18292,0.18833,0.1939,0.19964,0.20555,0.21163,0.21664,0.2213,0.22606,0.23092,0.23589,0.24096,0.24614,0.25009,0.25179,0.2535,0.25523,0.25697,0.25872,0.26048,0.26225,0.26404,0.26583,0.26764,0.26915,0.26913,0.26912,0.2691,0.26908,0.26906,0.26905,0.26903,0.26901,0.26899,0.26898,0.26896,0.26894,0.26892,0.26891,0.26888,0.26886,0.26883,0.26881,0.26878,0.26875,0.26873,0.2687,0.26868,0.26865,0.26862,0.2686,0.26857,0.26854,0.26747,0.2657,0.26394,0.26219,0.26045,0.25873,0.25702,0.25532,0.25258,0.24724,0.24202,0.2369,0.2319,0.22699,0.2222,0.2175,0.2129,0.20766,0.20199,0.1967,0.19155,0.18653,0.18165,0.17689,0.17226,0.16775,0.16336,0.15908,0.15491,0.15086,0.14691,0.14306,0.13931,0.13567,0.13211,0.12865,0.12528,0.122,0.11881,0.1157,0.11229,0.10863,0.10509,0.10166,0.098342,0.095133,0.092029,0.089027,0.086122,0.083313,0.080594,0.077965,0.075421,0.07526,0.076163,0.077078,0.078003,0.078939,0.079982,0.081135,0.082305,0.083492,0.084696,0.085918,0.087157,0.088414,0.089689,0.090982,0.092294,0.093625,0.094975,0.096345,0.097734,0.099143,0.10069,0.097109,0.093652,0.090318,0.087103,0.084003,0.081012,0.078129,0.075347,0.072665,0.070079,0.067584,0.065178,0.062696,0.060152,0.057711,0.055369,0.053122,0.050967,0.048899,0.046914,0.045011,0.043184,0.041432,0.039751,0.038039,0.036307,0.034654,0.033076,0.03157,0.030133,0.028761,0.027451,0.026201,0.025008,0.023869,0.022783,0.021745,0.020701,0.019657]';
            end

        case 13
            %just the modified middle ear transfer function from Glasberg
            %and Moore (2006) derived as above
            % values refined by Doheon Lee
            if faster == 0
                I=I.*[0.0024443,0.011856,0.023323,0.036415,0.045062,0.055763,0.069005,0.085392,0.10567,0.11946,0.13251,0.14699,0.16259,0.17354,0.18523,0.1977,0.21102,0.22523,0.23939,0.25333,0.26809,0.28371,0.30024,0.31773,0.33623,0.35421,0.36909,0.38459,0.40075,0.41758,0.42443,0.42984,0.43531,0.44086,0.44648,0.45217,0.45793,0.46376,0.46967,0.47566,0.48172,0.48786,0.49407,0.50037,0.5061,0.51161,0.51718,0.52281,0.5285,0.53425,0.54006,0.54594,0.55188,0.55788,0.56395,0.5663,0.56331,0.56033,0.55737,0.55443,0.5515,0.54859,0.54569,0.54281,0.53063,0.51522,0.50027,0.48575,0.47165,0.45796,0.44466,0.43176,0.41922,0.40705,0.39658,0.3878,0.37922,0.37082,0.36261,0.35459,0.34674,0.33906,0.33155,0.32422,0.31704,0.31002,0.30497,0.30007,0.29525,0.29051,0.28584,0.28125,0.27673,0.27229,0.26792,0.26361,0.25925,0.25491,0.25064,0.24644,0.24231,0.23825,0.23426,0.23034,0.22648,0.22269,0.21896,0.21529,0.21169,0.20814,0.20465,0.20123,0.19786,0.19454,0.19128,0.18808,0.18493,0.18183,0.17879,0.17579,0.17285,0.16995,0.16711,0.16431,0.16156,0.15885,0.15619,0.15357,0.151,0.14847,0.14599,0.14354,0.14114,0.13877,0.13645,0.13416,0.13192,0.12971,0.12753,0.1254,0.1233,0.12123,0.1192,0.11721,0.11524,0.11331,0.11141,0.10955,0.10771,0.10591,0.10414,0.10239,0.10068,0.09899,0.097333,0.095702,0.094099,0.095128,0.0964,0.097689,0.098995,0.10032,0.10166,0.10302,0.1044,0.10579,0.10721,0.10864,0.1101,0.11157,0.11306,0.11457,0.1161,0.11766,0.11923,0.12082,0.12244,0.12408,0.12574,0.12742,0.12912,0.13085,0.13675,0.1431,0.14974,0.1567,0.16397,0.17159,0.17955,0.18789,0.19143,0.19212,0.19281,0.1935,0.1942,0.1949,0.1956,0.19631,0.19702,0.19773,0.19844,0.19915,0.19987,0.20059,0.20131,0.20204,0.20277,0.2035,0.20423,0.20497,0.20571,0.20645,0.20719,0.20794,0.20869,0.20944,0.2102,0.21095,0.21171,0.21248,0.21324,0.21401,0.21478,0.21556,0.21633,0.21711,0.2179,0.21868,0.21947,0.22026,0.22105,0.22185,0.22225,0.22205,0.22186,0.22166,0.22146,0.22127,0.22107,0.22087,0.22068,0.22048,0.22029,0.22009,0.2199,0.2197,0.21951,0.21931,0.21912,0.21892,0.21873,0.21853,0.21834,0.21814,0.21795,0.21776,0.21756,0.21737,0.21718,0.21698,0.21679,0.2166,0.21641,0.21621,0.21602,0.21583,0.21564,0.21545,0.21526,0.21507,0.21487,0.21468,0.21449,0.2143,0.21411,0.21392,0.21373,0.21354,0.21335,0.21316,0.21297,0.21279,0.2126,0.21241,0.21187,0.21107,0.21027,0.20948,0.20868,0.20789,0.20711,0.20632,0.20554,0.20476,0.20399,0.20322,0.20245,0.20168,0.20092,0.20016,0.1994,0.19865,0.19789,0.19714,0.1964,0.19566,0.19491,0.19418,0.19344,0.19271,0.19198,0.19125,0.19053,0.18981,0.18909,0.18837,0.18766,0.18695,0.18624,0.18554,0.18484,0.18414,0.18344,0.18275,0.18205,0.18136,0.18068,0.17999,0.17931,0.17863,0.17796,0.17712,0.17561,0.17412,0.17264,0.17117,0.16972,0.16828,0.16685,0.16543,0.16402,0.16263,0.16125,0.15988,0.15852,0.15717,0.15583,0.15451,0.1532,0.1519,0.1506,0.14932,0.14806,0.1468,0.14555,0.14431,0.14309,0.14187,0.14066,0.13947,0.13828,0.13711,0.13594,0.13479,0.13364,0.13251,0.13138,0.13027,0.12916,0.12806,0.12697,0.12589,0.12482,0.12376,0.12271,0.12167,0.12063,0.11961,0.11859,0.11758,0.11659,0.11559,0.11461,0.11364,0.11267,0.11172,0.11077,0.10984,0.10891,0.108,0.10709,0.10619,0.10529,0.10441,0.10353,0.10266,0.1018,0.10094,0.10009,0.09925,0.098416,0.097588,0.096767,0.095954,0.095147,0.094346,0.093553,0.092766,0.091986,0.091213,0.090445,0.089685,0.088931,0.088183,0.087441,0.086706,0.085977,0.085254,0.084537,0.083826,0.083121,0.082422,0.081729,0.081041,0.08036,0.079684,0.079014,0.078349,0.07769,0.077037,0.076389,0.075747,0.07511,0.074478,0.073852,0.073231,0.072615,0.072004,0.071399,0.070798,0.070203,0.069612,0.069027,0.068446,0.067871,0.0673,0.066734,0.066303,0.066107,0.065911,0.065716,0.065521,0.065327,0.065133,0.06494,0.064748,0.064556,0.064364,0.064174,0.063983,0.063794,0.063605,0.063416,0.063228,0.063041,0.062854,0.062668,0.062482,0.062297,0.062112,0.061928,0.061745,0.061562,0.061379,0.061197,0.061016,0.060835,0.060655,0.060475,0.060296,0.060117,0.059939,0.059761,0.059584,0.059407,0.059231,0.059056,0.058881,0.058706,0.058532,0.058359,0.058186,0.058013,0.057841,0.05767,0.057499,0.057329,0.057159,0.057066,0.057373,0.057683,0.057994,0.058307,0.058622,0.058938,0.059256,0.059576,0.059897,0.060221,0.060546,0.060872,0.061201,0.061531,0.061863,0.062197,0.062533,0.06287,0.063209,0.06355,0.063893,0.064238,0.064585,0.064933,0.065284,0.065636,0.06599,0.066346,0.066704,0.067064,0.067426,0.06779,0.068156,0.068523,0.068893,0.069265,0.069639,0.070014,0.070392,0.070772,0.071154,0.071538,0.071924,0.072312,0.072702,0.073094,0.073489,0.073885,0.074284,0.074685,0.075088,0.075493,0.075901,0.07631,0.076722,0.077136,0.077552,0.077971,0.078391,0.078814,0.07924,0.079667,0.080097,0.080529,0.080964,0.081401,0.08184,0.082282,0.082726,0.083172,0.083621,0.084072,0.084526,0.084982,0.08544,0.085901,0.086124,0.086347,0.086571,0.086795,0.08702,0.087246,0.087472,0.087699,0.087926,0.088154,0.088383,0.088612,0.088841,0.089072,0.089302,0.089534,0.089766,0.089999,0.090232,0.090466,0.0907,0.090935,0.091171,0.091407,0.091644,0.091882,0.09212,0.092359,0.092598,0.092838,0.093079,0.09332,0.093562,0.093805,0.094048,0.094292,0.094536,0.094781,0.095027,0.095273,0.09552,0.095767,0.096016,0.096265,0.096514,0.096764,0.097015,0.097267,0.097519,0.097771,0.098025,0.098279,0.098032,0.097325,0.096624,0.095928,0.095236,0.09455,0.093869,0.093192,0.092521,0.091854,0.091192,0.090535,0.089883,0.089235,0.088592,0.087953,0.08732,0.08669,0.086066,0.085445,0.08483,0.084218,0.083611,0.083009,0.082411,0.081817,0.081227,0.080642,0.080061,0.079484,0.078911,0.078343,0.077778,0.077217,0.076661,0.076109,0.07556,0.075016,0.074475,0.073938,0.073406,0.072877,0.072351,0.07183,0.071312,0.070799,0.070288,0.069782,0.069279,0.06878,0.068284,0.067792,0.067303,0.066818,0.066337,0.065859,0.065384,0.064913,0.064445,0.063981,0.06352,0.063062,0.062608,0.062156,0.061708,0.061263,0.060822,0.060383,0.059948,0.059516,0.059087,0.058661,0.058238,0.057818,0.057401,0.056987,0.056576,0.056169,0.055764,0.055362,0.054962,0.054566,0.054173,0.053782,0.053395,0.05301,0.052627,0.052248,0.051871,0.051497,0.051126,0.050757,0.050392,0.050028,0.049668,0.049309,0.048954,0.048601,0.048251,0.047903,0.047557,0.047215,0.046874,0.046536,0.046201,0.045868,0.045537,0.045209,0.044883,0.044559,0.044238,0.043919,0.043602,0.043288,0.042976,0.042666,0.042359,0.042053,0.04175,0.041449,0.041152,0.040864,0.040579,0.040295,0.040013,0.039733,0.039456,0.03918,0.038906,0.038634,0.038364,0.038096,0.037829,0.037565,0.037302,0.037041,0.036782,0.036525,0.03627,0.036016,0.035765,0.035515,0.035266,0.03502,0.034775,0.034532,0.034342,0.034188,0.034035,0.033882,0.03373,0.033579,0.033428,0.033278,0.033129,0.03298,0.032832,0.032685,0.032538,0.032392,0.032247,0.032102,0.031958,0.031814,0.031672,0.03153,0.031388,0.031247,0.031107,0.030968,0.030829,0.03069,0.030537,0.030384,0.030233,0.030082,0.029932,0.029782,0.029633,0.029485,0.029338,0.029192,0.029046,0.028901,0.028757,0.028613,0.02847,0.028328,0.028186,0.028046,0.027906,0.027766,0.027628,0.02749,0.027352,0.027216,0.02708,0.026945,0.02681,0.026676,0.026543,0.026411,0.026279,0.026147,0.026017,0.025887,0.025758,0.025629,0.025501,0.025374,0.025247,0.025121,0.024995,0.024871,0.024746,0.024623,0.0245,0.024378,0.024256,0.024135,0.024014,0.023894,0.023775,0.023656,0.023538,0.023421,0.023304,0.023187,0.023072,0.022956,0.022842,0.022728,0.022614,0.022501,0.022389,0.022277,0.022166,0.022055,0.021945,0.021835,0.021726,0.021618,0.02151,0.021402,0.021296,0.021189,0.021083,0.020978,0.020873,0.020769,0.020665,0.020562,0.02046,0.020357,0.020256,0.020155,0.020054,0.019954,0.019854,0.019755,0.019656,0.019558,0.019461,0.019363,0.019267,0.019171,0.019075,0.01898,0.018885,0.01879,0.018697,0.018603,0.01851,0.018418,0.018326,0.018234,0.018143,0.018053,0.017963,0.017873,0.017784,0.017695,0.017607,0.017519,0.017431,0.017344,0.017257,0.017171,0.017086,0.017,0.016915,0.016831,0.016747,0.016663,0.01658,0.016497,0.016415,0.016333,0.016251,0.01617,0.016089,0.016009,0.015929,0.01585,0.01577,0.015692,0.015613,0.015535,0.015458,0.015381,0.015304,0.015227,0.015151,0.015076,0.015,0.014925,0.014851,0.014777,0.014703,0.01463,0.014557,0.014484,0.014412,0.01434,0.014268,0.014197,0.014126,0.014055,0.013985,0.013915,0.013846,0.013777]';

            elseif faster ==1
                I=I.*[0.0024443,0.011856,0.023323,0.036415,0.045062,0.055763,0.069005,0.085392,0.10567,0.11946,0.13251,0.14699,0.16259,0.17354,0.18523,0.1977,0.21102,0.22523,0.23939,0.25333,0.26809,0.28371,0.30024,0.31773,0.33623,0.35421,0.36909,0.38459,0.40075,0.41758,0.42443,0.42984,0.43531,0.44086,0.44648,0.45217,0.45793,0.46376,0.46967,0.47566,0.48172,0.48786,0.49407,0.50037,0.5061,0.51161,0.51718,0.52281,0.5285,0.53425,0.54006,0.54594,0.55188,0.55788,0.56395,0.5663,0.56331,0.56033,0.55737,0.55443,0.5515,0.54859,0.54569,0.54281,0.53063,0.51522,0.50027,0.48575,0.47165,0.45796,0.44466,0.43176,0.41922,0.40705,0.39658,0.3878,0.37922,0.37082,0.36261,0.35459,0.34674,0.33906,0.33155,0.32422,0.31704,0.31002,0.30497,0.30007,0.29525,0.29051,0.28584,0.28125,0.27673,0.27229,0.26792,0.26361,0.25925,0.25491,0.25064,0.24644,0.24231,0.23825,0.23426,0.23034,0.22648,0.22269,0.21896,0.21529,0.21169,0.20814,0.20465,0.20123,0.19786,0.19454,0.19128,0.18808,0.18493,0.18183,0.17879,0.17579,0.17285,0.16995,0.16711,0.16431,0.16156,0.15885,0.15619,0.15357,0.151,0.14847,0.14599,0.14354,0.14114,0.13877,0.13645,0.13416,0.13192,0.12971,0.12753,0.1254,0.1233,0.12123,0.1192,0.11721,0.11524,0.11331,0.11141,0.10955,0.10771,0.10591,0.10414,0.10239,0.10068,0.09899,0.097333,0.095702,0.094099,0.095128,0.0964,0.097689,0.098995,0.10032,0.10166,0.10302,0.1044,0.10579,0.10721,0.10864,0.1101,0.11157,0.11381,0.11688,0.12002,0.12326,0.12657,0.12998,0.13989,0.15318,0.16774,0.18367,0.19177,0.19316,0.19455,0.19595,0.19737,0.1988,0.20023,0.20168,0.20313,0.2046,0.20608,0.20757,0.20906,0.21057,0.2121,0.21363,0.21517,0.21672,0.21829,0.21987,0.22145,0.22215,0.22176,0.22137,0.22097,0.22058,0.22019,0.2198,0.21941,0.21902,0.21863,0.21824,0.21785,0.21747,0.21708,0.2167,0.21631,0.21583,0.21526,0.21468,0.21411,0.21354,0.21297,0.21241,0.21027,0.20789,0.20554,0.20322,0.20092,0.19865,0.1964,0.19418,0.19198,0.18981,0.18766,0.18554,0.18344,0.18136,0.17931,0.17712,0.17264,0.16828,0.16402,0.15988,0.15583,0.1519,0.14806,0.14431,0.14007,0.13537,0.13082,0.12643,0.12219,0.11809,0.11412,0.1103,0.10664,0.10309,0.09967,0.09636,0.093159,0.090064,0.087073,0.08418,0.081384,0.078681,0.076067,0.073541,0.071098,0.068736,0.066453,0.065521,0.064556,0.063605,0.062668,0.061745,0.060835,0.059939,0.059056,0.058186,0.057329,0.057994,0.059576,0.061201,0.06287,0.064585,0.066346,0.068156,0.070014,0.072118,0.074484,0.076929,0.079453,0.08206,0.084753,0.086683,0.08804,0.089418,0.090818,0.092239,0.093683,0.09515,0.096639,0.098152,0.094893,0.090535,0.086066,0.081817,0.077778,0.073938,0.070288,0.066818,0.06352,0.060383,0.057401,0.054566,0.051871,0.049309,0.046705,0.044078,0.041599,0.039318,0.037172,0.035143,0.033654,0.032465,0.031318,0.030157,0.028973,0.027836,0.026676,0.025501,0.024378,0.023304,0.022277,0.021296,0.020357,0.019461,0.018603,0.017784,0.017,0.016251,0.015535,0.014814,0.01409]';
            end

            % case 14
            % ADD MORE FILTERS HERE AS DESIRED
            % frequencies are:
            % ((3:961)-1)*32000/2048; and
            % [31.25:15.625:2671.875,2695.3125:31.25:4132.8125,4171.875:46.875:5578.125,5632.8125:62.5:7007.8125,7078.125:78.125:8406.25,8492.1875:93.75:9898.4375,10000:109.375:11312.5,11429.6875:125:12804.6875,12937.5:140.625:14625,14773.4375,14929.6875];
    end %switch filtermethod

    % use if you wish to monitor the sound pressure level here
    % level = 10*log10(sum(I)+eps)
    for channel=1:nchan
        % roex component levels
        % this is slow due to the large number of spectrum components
        % tic
        intensity = zeros(nfreq,nfreq);
        intensity(indg) = (1+p(col)'.*(g(indg))).* exp(-p(col)'.*(g(indg))) .* I(row,channel);
        RoexL = 10.*log10(sum(intensity,1)+eps);
        % toc

        % excitation pattern
        p2(g2neg) = p51(g2neg)-0.35.*(p51(g2neg)./p511k) .* (RoexL(row_g2_neg(:))-51)';
        if p2(g2lessthan2) <0.1, p2(g2lessthan2) = 0.1; end
        Imat = repmat(I(:,channel),1,nERBS);
        Esig = sum(g2lessthan2.*Imat .* (1+p2.*abs(g2)).*exp(-p2.*abs(g2)))'+eps;
        Elevel = 10*log10(Esig);

        % loudness
        if nchan == 1
            Highlevels=Elevel>100;
            SpecLoud(Highlevels) = C .* (Esig(Highlevels) ./ 1.04e6).^0.5;
            Midlevels = false(nERBS,1);
            Midlevels(~Highlevels)=Esig(~Highlevels)>Ethrq(~Highlevels);
            SpecLoud(Midlevels) = C .* (((g_N(Midlevels) .* Esig(Midlevels)+ A(Midlevels)).^alpha(Midlevels))...
                - (A(Midlevels).^alpha(Midlevels)));
            Lowlevels = ~(Highlevels | Midlevels);
            SpecLoud(Lowlevels)= C .* (((2 .* Esig(Lowlevels)) ./ (Esig(Lowlevels) + Ethrq(Lowlevels))).^1.5)...
                .* (((g_N(Lowlevels) .* Esig(Lowlevels) + A(Lowlevels)).^alpha(Lowlevels)) - (A(Lowlevels).^alpha(Lowlevels)));
            SpecLoud = SpecLoud .* 2;
        else
            if channel == 1
                Highlevels=Elevel>100;
                SpecLoudL(Highlevels) = C .* (Esig(Highlevels) ./ 1.04e6).^0.5;
                Midlevels = false(nERBS,1);
                Midlevels(~Highlevels)=Esig(~Highlevels)>Ethrq(~Highlevels);
                SpecLoudL(Midlevels) = C .* (((g_N(Midlevels) .* Esig(Midlevels)+ A(Midlevels)).^alpha(Midlevels))...
                    - (A(Midlevels).^alpha(Midlevels)));
                Lowlevels = ~(Highlevels | Midlevels);
                SpecLoudL(Lowlevels)= C .* (((2 .* Esig(Lowlevels)) ./ (Esig(Lowlevels) + Ethrq(Lowlevels))).^1.5)...
                    .* (((g_N(Lowlevels) .* Esig(Lowlevels) + A(Lowlevels)).^alpha(Lowlevels)) - (A(Lowlevels).^alpha(Lowlevels)));

                SLsmooth=repmat(SpecLoudL, 1, nERBS).*w_bi;
                SLsmooth=sum(SLsmooth, 2);
            elseif channel == 2
                Highlevels=Elevel>100;
                SpecLoudR(Highlevels) = C .* (Esig(Highlevels) ./ 1.04e6).^0.5;
                Midlevels = false(nERBS,1);
                Midlevels(~Highlevels)=Esig(~Highlevels)>Ethrq(~Highlevels);
                SpecLoudR(Midlevels) = C .* (((g_N(Midlevels) .* Esig(Midlevels)+ A(Midlevels)).^alpha(Midlevels))...
                    - (A(Midlevels).^alpha(Midlevels)));
                Lowlevels = ~(Highlevels | Midlevels);
                SpecLoudR(Lowlevels)= C .* (((2 .* Esig(Lowlevels)) ./ (Esig(Lowlevels) + Ethrq(Lowlevels))).^1.5)...
                    .* (((g_N(Lowlevels) .* Esig(Lowlevels) + A(Lowlevels)).^alpha(Lowlevels)) - (A(Lowlevels).^alpha(Lowlevels)));

                SRsmooth=repmat(SpecLoudR, 1, nERBS).*w_bi;
                SRsmooth=sum(SRsmooth,2);
                % inhibit specific loudness
                INHIBL=2./(1+(sech((SRsmooth+10e-13)./(SLsmooth+10e-13))).^p_bi);
                INHIBR=2./(1+(sech((SLsmooth+10e-13)./(SRsmooth+10e-13))).^p_bi);
                SpecLoud=SpecLoudL./INHIBL + SpecLoudR./INHIBR;
            end %if channel ==
        end % if nchan
    end % for channel = 1:nchan
    index = windowcount+1;
    InstantaneousLoudness(index) = sum(SpecLoud)*0.25; % multiply by erb step
    if InstantaneousLoudness(index) < minloud
        InstantaneousLoudness(index) = minloud;
    end

    %Short-term loudness
    if  InstantaneousLoudness(index) > ShortTermLoudness(index-1)
        ShortTermLoudness(index)=XSa*InstantaneousLoudness(index) +...
            (1-XSa)*ShortTermLoudness(index-1);
    else
        ShortTermLoudness(index)=XSr*InstantaneousLoudness(index) +...
            (1-XSr)*ShortTermLoudness(index-1);
    end %if

    %Long-term loudness
    if  ShortTermLoudness(index) > LongTermLoudness(index-1)
        LongTermLoudness(index)=XLa*ShortTermLoudness(index) +...
            (1-XLa)*LongTermLoudness(index-1);
    else
        LongTermLoudness(index)=XLr*ShortTermLoudness(index) +...
            (1-XLr)*LongTermLoudness(index-1);
    end %if

end % for windowcount = 1:nwindows

% clear some memory in case it is needed for the figure
clear intensity g g2 p2

% additional silence at the end, to calculate loudness decay
if decay > 0
    InstantaneousLoudness = [InstantaneousLoudness;zeros(decay,1).*minloud];
    ShortTermLoudness = [ShortTermLoudness;zeros(decay,1)];
    LongTermLoudness = [LongTermLoudness;zeros(decay,1)];
    for k=1:decay+1
        index = nwindows+1+k;
        %Short-term loudness
        if  InstantaneousLoudness(index) > ShortTermLoudness(index-1)
            ShortTermLoudness(index)=XSa*InstantaneousLoudness(index) +...
                (1-XSa)*ShortTermLoudness(index-1);
        else
            ShortTermLoudness(index)=XSr*InstantaneousLoudness(index) +...
                (1-XSr)*ShortTermLoudness(index-1);
        end %if

        %Long-term loudness
        if  ShortTermLoudness(index) > LongTermLoudness(index-1)
            LongTermLoudness(index)=XLa*ShortTermLoudness(index) +...
                (1-XLa)*LongTermLoudness(index-1);
        else
            LongTermLoudness(index)=XLr*ShortTermLoudness(index) +...
                (1-XLr)*LongTermLoudness(index-1);
        end %if
    end % for k=1:decay+1
end % if decay > 0

switch doplot
    case 1
        % plot of loudness
        figure
        plot(times,InstantaneousLoudness,'k');
        hold on
        plot(times,ShortTermLoudness,'r');
        plot(times,LongTermLoudness,'b');
        xlabel('Time (s)');
        ylabel('Loudness (sone)');
        hold off

    case 2 | 3
        % plot of loudness level
        % combine loudness into matrix for conversion to loudness level
        N=[InstantaneousLoudness, ShortTermLoudness, LongTermLoudness];

        % The following sone to phon conversion was derived from this loudness model
        % by calculating the loudness a 1 kHz tone at many sound pressure levels
        % (from -20 dB to 139 dB, at 0.1 dB steps). The tone was a single
        % spectrum component as it entered the loudness model (i.e. prior to
        % the roex calculation), and so was unaffected by the FFT
        % analysis and outer/middle ear stages of the model. The time-varying
        % and binaural stages of the model were not used either. The polynomials below
        % are based on these calculated values, and are accurate within a +/- 0.5
        % phon range.
        
        Nsmall = round(N.*1e4)/1e4 <= minloud; 
        Nlarge = N > 4.6346e+03; %values that are out of range (greater than 130 dB)
        Medium_crossover=4.181876451206376; % a crossover point for applying the two polynominal function
        Nmedium1 = minloud < round(N.*1e4)/1e4 & N < Medium_crossover;
        Nmedium2 = Medium_crossover <= N  & N < 4.6346e+03;

        log2sone = zeros(nwindows+2+decay,3);
        log2sone(Nmedium1) = log2(N(Nmedium1));
        log2sone(Nmedium2) = log2(N(Nmedium2));

        LN=zeros(nwindows+2+decay, 3); %loudness level
        LN(Nsmall) = 2;
        LN(Nlarge)=NaN;
        LN(Nmedium1) = -0.000001891383225 .* log2sone(Nmedium1).^6 -0.000131861311930 .* log2sone(Nmedium1).^5  ...
            -0.002885410980939 .* log2sone(Nmedium1).^4 -0.006304645634479 .* log2sone(Nmedium1).^3 ...
            +0.577697444612570 .* log2sone(Nmedium1).^2 +8.696620195281991 .* log2sone(Nmedium1) ...
            +40; %+39.914029611835041;
        LN(Nmedium2) =  -0.000271548662498 .* log2sone(Nmedium2).^6 +0.011957866530551 .* log2sone(Nmedium2).^5 ...
            -0.201870560271078 .* log2sone(Nmedium2).^4  +1.636339276141087 * log2sone(Nmedium2).^3 ...
            -6.997514386301266  * log2sone(Nmedium2).^2 +25.386669115791680 * log2sone(Nmedium2) ...
            +26.087617846072916; %+26.472581191744133;
        
        if doplot == 3
            InstantaneousLoudness = LN(:,1);
            ShortTermLoudness = LN(:,2);
            LongTermLoudness = LN(:,3);
        end % if doplot == 3

        clear log2sone
        figure
        plot(times,LN(:,1),'k')
        hold on
        plot(times,LN(:,2),'r')
        plot(times,LN(:,3),'b')
        xlabel('Time (s)')
        ylabel('Loudness Level (phon)')
        hold off
end; % switch doplot