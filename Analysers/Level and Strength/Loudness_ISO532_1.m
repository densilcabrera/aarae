function [OUT] = Loudness_ISO532_1(IN, field, method, TIME_SKIP)
% This function implements Fastl & Zwicker Loudness model for stationary 
% signals (Method A) and for arbitrary signals (Method B) as detailed in
% ISO 532-1 "Acousics - Methods for calculating loudness - Part 1: Zwicker
% method", 2014.
%
% INPUT ARGUMENTS
% IN - audio signal (1 channel only as specified by the standard)
% field - free field = 0; diffuse field = 1
% method - stationary = 0; time varying = 1;
% TIME_SKIP - in seconds for level calculation (stationary signals) 
%
% OUTPUTS
% * Time-varying Loudness
% * Time-averaged Specific Loudness
% * Time-averaged Loudness statistics
% Source: C code is provided in the ISO532 Annex A (2014)
% MATLAB code with adjustments to the C code by Ella Manor for AARAE (2015)

% *************************************************************************

if nargin == 1 

    param = inputdlg({'Sound field: Free (0); Diffuse (1)';... 
                      'Method: Stationary (0); Time varying (1)';...
                      'Start calculation at: (seconds)'},...
                      'Parameters',... 
                      [1 50],... 
                      {'0';'0';'0.'}); % And the preset answers for your dialog.
    if str2double(param(2)) == 1 && str2double(param(3)) > 0
        h=warndlg('When using time varying method the calculation must start from time 0');
        uiwait(h)
        param{3} = '0';
    end
    
    if length(param) < 2, param = []; end % You should check that the user 
                                          % has input all the required
                                          % fields.
    if ~isempty(param) % If they have, you can then assign the dialog's
                       % inputs to your function's input parameters.
        field = str2double(param(1));
        method = str2double(param(2));
        TIME_SKIP = str2double(param(3));
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
    
else
    param = [];
end

% *************************************************************************
if isstruct(IN) 
    
    audio = IN.audio; 
    fs = IN.fs; 
    
    if (fs ~= 32000)
        if (fs ~= 44100) 
            if (fs ~= 48000)
                h=warndlg('Sampling rates supported by this analyser include: 32 kHz, 44.1 kHz and 48 kHz');
                uiwait(h)
                OUT = IN;
                return % get out of here!
            end
        end
    end
    
    if isfield(IN,'cal') % Get the calibration offset if it exists
        cal = IN.cal;
    else
        h=warndlg('Calibration data missing - please calibrate now.','Name of Analyser','modal');
        uiwait(h)
        IN = cal_aarae(IN); % cal_aarae is a Basic Processor that emulates AARAE's 'CAL' button
        if isempty(IN) % user pressed cancel within cal_aarae, or the calibration failed
            OUT = IN; % you need to return an empty output
            return % get out of here!
        end
        cal = IN.cal; % get the cal field values
    end
    
    if isfield(IN,'name') % Get the AARAE name if it exists
        name = IN.name; % this is a string
    else
        name = ''; % empty string - can be concatenated without any problems
    end
end

% *************************************************************************

if ~isempty(audio) && ~isempty(fs)
    
    calconstant = 0; % 1 kHz at 40 dB = 1 sone
    cal = cal-calconstant;
    if length(cal) == 1
        audio = audio .* 10.^(cal/20);
    else
        audio(:,1) = audio(:,1) .* 10.^(cal(1)/20);
        audio(:,2) = audio(:,2) .* 10.^(cal(2)/20);
    end

    disp(['rms level of the entire wave ', num2str(10*log10(mean(audio.^2)+10e-99)+calconstant), ' dB'])

    [len,chans] =  size(audio);   
    
    if chans > 1
        audio = audio(:,1);
        disp('Only one channel signals are analysed by Loudness_ISO532.') 
        chans = 1;
    end
    
% Time constants for non-linear temporal decay
Tshort = 0.005;
Tlong = 0.015;
Tvar = 0.075;

% Factors for virtual upsampling/inner iterations
NL_ITER = 24;
% Sampling rate to which third-octave-levels are downsampled 
SR_LEVEL = 2000;
% Sampling rate to which output/total summed loudness is downsampled 
SR_LOUDNESS = 500;  
% Tiny value for adjusting intensity levels for stationary signals
TINY_VALUE = 1e-12; 
% ref value for SPL calculation
I_REF = 4e-10; 

% ***************************
% STEP 1 - Resample to 48 kHz
% ***************************

if fs ~= 48000
    gcd_fs = gcd(48000,fs); % greatest common denominator
    audio = resample(audio,48000/gcd_fs,fs/gcd_fs);
    fs = 48000;
    len = size(audio,1);
end

% ******************************************************************
% Assign values to global variables according to the selected method
% ******************************************************************

if method == 0
    SampleRateLevel = 1;
    DecFactorLevel = len;
    DecFactorLoudness = 1;
    NumSamplesLevel = 1;
elseif method == 1
    SampleRateLevel = SR_LEVEL;
    SampleRateLoudness = SR_LOUDNESS;
    DecFactorLevel = fs/SampleRateLevel;
    DecFactorLoudness = SampleRateLevel/SampleRateLoudness;
    NumSamplesLevel = ceil(len/DecFactorLevel);
    NumSamplesLoudness = ceil(NumSamplesLevel/DecFactorLoudness);

end

% **************************************************
% STEP 2 - Create filter bank and filter the signal
% **************************************************

% reference 
br = [1,2,1;1,0,-1;1,-2,1];
ar = [1,-2,1;1,-2,1;1,-2,1];

CenterFrequency = zeros(28,1);
for i = 1:28
    CenterFrequency(i) = 10^(((i-1)-16)/10.) * 1000; % calculate centre frequencies
end

% filter 'a' coefficient offsets TABLES A.1 A.2
ad = cat(3,[0,-6.70260e-004,6.59453e-004;...
0,-3.75071e-004,3.61926e-004;...
0,-3.06523e-004,2.97634e-004],... % 25 Hz
[0,-8.47258e-004,8.30131e-004;...
0,-4.76448e-004,4.55616e-004;...
0,-3.88773e-004,3.74685e-004],... % 31.5 Hz
[0,-1.07210e-003,1.04496e-003;...
0,-6.06567e-004,5.73553e-004;...
0,-4.94004e-004,4.71677e-004],... % 40 Hz
[0,-1.35836e-003,1.31535e-003;...
0,-7.74327e-004,7.22007e-004;...
0,-6.29154e-004,5.93771e-004],... % 50 Hz
[0,-1.72380e-003,1.65564e-003;...
0,-9.91780e-004,9.08866e-004;...
0,-8.03529e-004,7.47455e-004],... % 63 Hz
[0,-2.19188e-003,2.08388e-003;...
0,-1.27545e-003,1.14406e-003;...
0,-1.02976e-003,9.40900e-004],... % 80 Hz
[0,-2.79386e-003,2.62274e-003;...
0,-1.64828e-003,1.44006e-003;...
0,-1.32520e-003,1.18438e-003],... % 100 Hz
[0,-3.57182e-003,3.30071e-003;...
0,-2.14252e-003,1.81258e-003;...
0,-1.71397e-003,1.49082e-003],... % 125 Hz
[0,-4.58305e-003,4.15355e-003;...
0,-2.80413e-003,2.28135e-003;...
0,-2.23006e-003,1.87646e-003],... % 160 Hz
[0,-5.90655e-003,5.22622e-003;...
0,-3.69947e-003,2.87118e-003;...
0,-2.92205e-003,2.36178e-003],... % 200 Hz
[0,-7.65243e-003,6.57493e-003;...
0,-4.92540e-003,3.61318e-003;...
0,-3.86007e-003,2.97240e-003],... % 250 Hz
[0,-1.00023e-002,8.29610e-003;...
0,-6.63788e-003,4.55999e-003;...
0,-5.15982e-003,3.75306e-003],... % 315 Hz
[0,-1.31230e-002,1.04220e-002;...
0,-9.02274e-003,5.73132e-003;...
0,-6.94543e-003,4.71734e-003],... % 400 Hz
[0,-1.73693e-002,1.30947e-002;...
0,-1.24176e-002,7.20526e-003;...
0,-9.46002e-003,5.93145e-003],... % 500 Hz
[0,-2.31934e-002,1.64308e-002;...
0,-1.73009e-002,9.04761e-003;...
0,-1.30358e-002,7.44926e-003],... % 630 Hz
[0,-3.13292e-002,2.06370e-002;...
0,-2.44342e-002,1.13731e-002;...
0,-1.82108e-002,9.36778e-003],... % 800 Hz
[0,-4.28261e-002,2.59325e-002;...
0,-3.49619e-002,1.43046e-002;...
0,-2.57855e-002,1.17912e-002],... % 1000 Hz
[0,-5.91733e-002,3.25054e-002;...
0,-5.06072e-002,1.79513e-002;...
0,-3.69401e-002,1.48094e-002],... % 1250 Hz
[0,-8.26348e-002,4.05894e-002;...
0,-7.40348e-002,2.24476e-002;...
0,-5.34977e-002,1.85371e-002],... % 1600 Hz
[0,-1.17018e-001,5.08116e-002;...
0,-1.09516e-001,2.81387e-002;...
0,-7.85097e-002,2.32872e-002],... % 2000 Hz
[0,-1.67714e-001,6.37872e-002;...
0,-1.63378e-001,3.53729e-002;...
0,-1.16419e-001,2.93723e-002],... % 2500 Hz
[0,-2.42528e-001,7.98576e-002;...
0,-2.45161e-001,4.43370e-002;...
0,-1.73972e-001,3.70015e-002],... % 3150 Hz
[0,-3.53142e-001,9.96330e-002;...
0,-3.69163e-001,5.53535e-002;...
0,-2.61399e-001,4.65428e-002],... % 4000 Hz
[0,-5.16316e-001,1.24177e-001;...
0,-5.55473e-001,6.89403e-002;...
0,-3.93998e-001,5.86715e-002],... % 5000 Hz
[0,-7.56635e-001,1.55023e-001;...
0,-8.34281e-001,8.58123e-002;...
0,-5.94547e-001,7.43960e-002],... % 6300 Hz
[0,-1.10165e+000,1.91713e-001;...
0,-1.23939e+000,1.05243e-001;...
0,-8.91666e-001,9.40354e-002],... % 8000 Hz
[0,-1.58477e+000,2.39049e-001;...
0,-1.80505e+000,1.28794e-001;...
0,-1.32500e+000,1.21333e-001],... % 10000 Hz
[0,-2.50630e+000,1.42308e-001;...
0,-2.19464e+000,2.76470e-001;...
0,-1.90231e+000,1.47304e-001]); % 12500 Hz

% filter gains
filtgain = [4.30764e-011;... % 25 Hz
8.59340e-011;... % 31.5 Hz
1.71424e-010;... % 40 Hz
3.41944e-010;... % 50 Hz
6.82035e-010;... % 63 Hz
1.36026e-009;... % 80 Hz
2.71261e-009;... % 100 Hz
5.40870e-009;... % 125 Hz
1.07826e-008;... % 160 Hz
2.14910e-008;... % 200 Hz
4.28228e-008;... % 250 Hz
8.54316e-008;... % 315 Hz
1.70009e-007;... % 400 Hz
3.38215e-007;... % 500 Hz
6.71990e-007;... % 630 Hz
1.33531e-006;... % 800 Hz
2.65172e-006;... % 1000 Hz
5.25477e-006;... % 1250 Hz
1.03780e-005;... % 1600 Hz
2.04870e-005;... % 2000 Hz
4.05198e-005;... % 2500 Hz
7.97914e-005;... % 3150 Hz
1.56511e-004;... % 4000 Hz
3.04954e-004;... % 5000 Hz
5.99157e-004;... % 6300 Hz
1.16544e-003;... % 8000 Hz
2.27488e-003;... % 10000 Hz
3.91006e-003];   % 12500 Hz

filteredaudio = zeros(len,28);


for n = 1:28
    filteredaudio(:,n) = filtgain(n) * filter(br(3,:),ar(3,:)-ad(3,:,n),...
    filter(br(2,:),ar(2,:)-ad(2,:,n),...
    filter(br(1,:),ar(1,:)-ad(1,:,n),audio)));
end

% **************************************************************
% STEP 3 - Squaring and smoothing by 3 1st order lowpass filters 
% **************************************************************

% square the filtered audio for both methods
filteredaudio = filteredaudio.^2;

NumSkip = floor(TIME_SKIP * fs);
smoothedaudio = zeros(len-NumSkip,28);
ThirdOctaveLevel = zeros(NumSamplesLevel,28);

if method == 1
    for i = 1:28
            if CenterFrequency(i) <= 1000
                Tau = 2/(3*CenterFrequency(i));
            else
                Tau = 2/(3*1000.);
            end

            % 3x smoothing 1st order low-pass filters in series
            A1 = exp(-1 ./ (fs * Tau));
            B0 = 1 - A1;
        smoothedaudio(:,i)= filter(B0,A1,...
            filter(B0,A1,...
            filter(B0,A1,filteredaudio(:,i))));
    end

    for i = 1:28
        soomthedaudioDecimated = decimate(filteredaudio(:,i),DecFactorLevel/6,12);
        soomthedaudioDecimated = decimate(soomthedaudioDecimated,DecFactorLevel/4,12);
        ThirdOctaveLevel(:,i) = real(10*log10(soomthedaudioDecimated+TINY_VALUE/I_REF));
    end
          
elseif method == 0
        if NumSkip > len/2
            warndlg('time signal too short');
            uiwait(h)
        end
        if NumSkip == 0; NumSkip = 1;end
        for i = 1:28
            smoothedaudio(1:len-NumSkip,i) = filteredaudio(NumSkip:len-1,i);
            ThirdOctaveLevel(NumSamplesLevel,i) = 10*log10((sum(smoothedaudio(:,i))/len)+TINY_VALUE/I_REF);
        end
end

% ***********************************************************
% STEP 4 - Apply weighting factor to first 3 1/1 octave bands
% ***********************************************************

% WEIGHTING BELLOW 315Hz TABLE A.3

% Ranges of 1/3 Oct bands for correction at low frequencies according to equal loudness contours
RAP = [45 55 65 71 80 90 100 120];

% Reduction of 1/3 Oct Band levels at low frequencies according to equal loudness contours 
% within the eight ranges defined by RAP (DLL)
DLL = [-32 -24 -16 -10 -5 0 -7 -3 0 -2 0;
    -29 -22 -15 -10 -4 0 -7 -2 0 -2 0;
    -27 -19 -14 -9 -4 0 -6 -2 0 -2 0;
    -25 -17 -12 -9 -3 0 -5 -2 0 -2 0;
    -23 -16 -11 -7 -3 0 -4 -1 0 -1 0;
    -20 -14 -10 -6 -3 0 -4 -1 0 -1 0;
    -18 -12 -9 -6 -2 0 -3 -1 0 -1 0;
    -15 -10 -8 -4 -2 0 -3 -1 0 -1 0];

CorrLevel = zeros(NumSamplesLevel,11);
Intens = zeros(NumSamplesLevel,11);
CBI = zeros(NumSamplesLevel,3);
LCB = zeros(NumSamplesLevel,3);

for j = 1:NumSamplesLevel
    for i = 1:11
        k=1;
        if (ThirdOctaveLevel(j,i) > (RAP(k)-DLL(k,i))) && (k < 8);
            k=k+1;
        end
        CorrLevel(j,i) = ThirdOctaveLevel(j,i) + DLL(k,i); % attenuated levels 
        Intens(j,i) = 10^(CorrLevel(j,i)/10); % attenuated 1/3 octave intensities
    end
    
    % *************************************************************
    % STEP 5 - Sumup intensity values of the first 3 critical bands
    % *************************************************************

    CBI(j,1) = sum(Intens(j,1:6)); % first critical band (sum of octaves (25Hz to 80Hz))
    CBI(j,2) = sum(Intens(j,7:9)); % second critical band (sum of octaves (100Hz to 160Hz))
    CBI(j,3) = sum(Intens(j,10:11)); % third critical band (sum of octaves (200Hz to 250Hz))
    
    FNGi = 10*log10(CBI);
    
    for i = 1:3
        if CBI(j,i)>0;
            LCB(j,i) = FNGi(j,i);
        else
            LCB(j,i) = 0;
        end
    end
end


% **********************************************************************
% STEP 6 - Calculate core loudness for each critical band
% **********************************************************************

% LEVEL CORRECTIONS TABLE A.5 (LDF0) DDF
% Level correction to convert from a free field to a diffuse field (last critical band 12.5kHz is not included)
DDF = [0 0 .5 .9 1.2 1.6 2.3 2.8 3 2 0 -1.4 -2 -1.9 -1 .5 3 4 4.3 4];

% LEVEL CORRECTIONS TABLE A.6 (LTQ)
% Critical band level at absolute threshold without taking into account the 
% transmission characteristics of the ear
LTQ = [30 18 12 8 7 6 5 4 3 3 3 3 3 3 3 3 3 3 3 3]; % Threshold due to internal noise
% Hearing thresholds for the excitation levels (each number corresponds to a critical band 12.5kHz is not included)

% LEVEL CORRECTIONS TABLE A.7 (LCB) DCB
% Correction factor because using third octave band levels (rather than critical bands)
DCB = [-.25 -.6 -.8 -.8 -.5 0 .5 1.1 1.5 1.7 1.8 1.8 1.7 1.6 1.4 1.2 .8 .5 0 -.5];

% LEVEL CORRECTIONS TABLE A.4 (A0)
% % Attenuation due to transmission in the middle ear
A0  = [0 0 0 0 0 0 0 0 0 0 -.5 -1.6 -3.2 -5.4 -5.6 -4 -1.5 2 5 12]; 
% Moore et al disagrees with this being flat for low frequencies

Le = zeros(NumSamplesLevel,20);
CoreL = zeros(NumSamplesLevel,21);

for j = 1:NumSamplesLevel
    for i = 1:20
        Le(j,i) = ThirdOctaveLevel(j,i+8);
        if i <= 3
            Le(j,i) = LCB(j,i);
        end
        Le(j,i) = Le(j,i) - A0(i);
        if field == 1;
            Le(j,i) = Le(j,i) + DDF(i);
        end
        if Le(j,i) > LTQ(i)
            S = 0.25;
            Le(j,i) = Le(j,i) - DCB(i);
            MP1 = 0.0635 * 10.^(0.025 * LTQ(i));
            MP2 = (1. - S + S*10^(0.1*(Le(j,i)-LTQ(i))))^0.25 - 1.;
            CoreL(j,i) = MP1 * MP2;
            if CoreL(j,i) <= 0
                CoreL(j,i) = 0;
            end
        end
    end
end


% *************************************************************************
% STEP 7 - Correction of specific loudness within the lowest critical band
% *************************************************************************

for j = 1:NumSamplesLevel
        CorrCL = 0.4 + 0.32 * CoreL(j,1)^.2;
        if CorrCL > 1; CorrCL = 1; end 
        CoreL(j,1) = CoreL(j,1)*CorrCL; 
end


% **********************************************************************
% STEP 8 - Implementation of NL Block
% **********************************************************************

if method == 1

    DeltaT = 1 / fs*NL_ITER;
    P = (Tvar + Tlong) / (Tvar*Tshort);
    Q = 1/(Tshort*Tvar);
    Lambda1 =-P/2 + sqrt(P*P/4 - Q);
    Lambda2 =-P/2 - sqrt(P*P/4 - Q);
    Den = Tvar * (Lambda1 - Lambda2);
    E1 = exp(Lambda1 * DeltaT);
    E2 = exp(Lambda2 * DeltaT);
    NlLpB(1) = (E1 - E2) / Den;
    NlLpB(2) =((Tvar * Lambda2 + 1) * E1 - (Tvar * Lambda1 + 1) * E2) / Den;
    NlLpB(3) =((Tvar * Lambda1 + 1) * E1 - (Tvar * Lambda2 +1 ) * E2) / Den;
    NlLpB(4) = (Tvar * Lambda1+1) * (Tvar * Lambda2 + 1) * (E1-E2) / Den;
    NlLpB(5) = exp(-DeltaT / Tlong);
    NlLpB(6) = exp(-DeltaT / Tvar);
    
    NlLpUoLast = 0;
    NlLpU2Last = 0;
    
    for i = 1:21
        for j = 1:NumSamplesLevel-1
            NextInput = CoreL(j+1,i);
            % interpolation steps between current and next sample
            Delta = (NextInput - CoreL(j,i)) / NL_ITER;
            Ui = CoreL(j,i);
            
            % f_nl_lp FUNCTION STARTS
            % case 1
            if Ui < NlLpUoLast
                if NlLpUoLast > NlLpU2Last
                    % case 1.1
                    U2 = NlLpUoLast*NlLpB(1) - NlLpU2Last*NlLpB(2);
                    Uo = NlLpUoLast*NlLpB(3) - NlLpU2Last*NlLpB(4);
                    if Ui > Uo
                        Uo  = Ui;
                    end
                    if U2 > Uo
                        U2 = Uo;
                    end
                else
                    % case 1.2
                    Uo = NlLpUoLast*NlLpB(5);
                    if Ui > Uo
                        Uo = Ui;
                    end
                    U2 = Uo;
                end
            % case 2
            elseif Ui == NlLpUoLast
                Uo = Ui;
                % case 2.1
                if Uo > NlLpUoLast
                    U2 = (NlLpUoLast - Ui)*NlLpB(6) + Ui;
                    % case 2.2
                else
                    U2 = Ui;
                end
            % case 3
            else
                Uo = Ui;
                U2 = (NlLpU2Last - Ui)*NlLpB(6) + Ui;
            end
            
            NlLpUoLast = Uo;
            NlLpU2Last = U2;
            
            CoreL(j,i) = Uo;
            % f_nl_lp FUNCTION ENDS

            Ui = Ui + Delta;
        
            % inner iteration
            for k = 1:NL_ITER
                % f_nl_lp FUNCTION STARTS
                % case 1
                if Ui < NlLpUoLast
                    if NlLpUoLast > NlLpU2Last
                        % case 1.1
                        U2 = NlLpUoLast*NlLpB(1) - NlLpU2Last*NlLpB(2);
                        Uo = NlLpUoLast*NlLpB(3) - NlLpU2Last*NlLpB(4);
                        if Ui > Uo
                            Uo  = Ui;
                        end
                        if U2 > Uo
                            U2 = Uo;
                        end
                    else
                        % case 1.2
                        Uo = NlLpUoLast*NlLpB(5);
                        if Ui > Uo
                            Uo = Ui;
                        end
                        U2 = Uo;
                    end
                % case 2
                elseif Ui == NlLpUoLast
                    Uo = Ui;
                    % case 2.1
                    if Uo > NlLpUoLast
                        U2 = (NlLpUoLast - Ui)*NlLpB(6) + Ui;
                        % case 2.2
                    else
                        U2 = Ui;
                    end
                % case 3
                else
                    Uo = Ui;
                    U2 = (NlLpU2Last - Ui)*NlLpB(6) + Ui;
                end

                NlLpUoLast = Uo;
                NlLpU2Last = U2;
                
                CoreL(j,i) = Uo;
                % f_nl_lp FUNCTION ENDS 
                Ui = Ui + Delta;
            end
        end
%      *d = f_nl_lp(*pCoreL, &NlLp);   % implemented in the source code
    end

end


% **********************************************************************
% STEP 9 - CALCULATE THE SLOPES
% **********************************************************************

% Upper limits of the approximated critical bands in Bark
% TABLE A.8
ZUP  = [.9 1.8 2.8 3.5 4.4 5.4 6.6 7.9 9.2 10.6 12.3 13.8 15.2 16.7 18.1 19.3 20.6 21.8 22.7 23.6 24];

% TABLE A.9
% Range of specific loudness for the determination of the steepness of the upper slopes in the specific loudness 
% - critical band rate pattern (used to plot the correct USL curve)
RNS = [21.5 18 15.1 11.5 9 6.1 4.4 3.1 2.13 1.36 .82 .42 .30 .22 .15 .10 .035 0];
% This is used to design the right hand slope of the loudness
USL = [13 8.2 6.3 5.5 5.5 5.5 5.5 5.5;
   9   7.5 6   5.1 4.5 4.5 4.5 4.5;
   7.8 6.7 5.6 4.9 4.4 3.9 3.9 3.9;
   6.2 5.4 4.6 4.0 3.5 3.2 3.2 3.2;
   4.5 3.8 3.6 3.2 2.9 2.7 2.7 2.7;
   3.7 3.0 2.8 2.35 2.2 2.2 2.2 2.2;
   2.9 2.3 2.1 1.9 1.8 1.7 1.7 1.7;
   2.4 1.7 1.5 1.35 1.3 1.3 1.3 1.3;
   1.95 1.45 1.3 1.15 1.1 1.1 1.1 1.1;
   1.5 1.2 .94 .86 .82 .82 .82 .82;
   .72 .67 .64 .63 .62 .62 .62 .62;
   .59 .53 .51 .50 .42 .42 .42 .42;
   .40 .33 .26 .24 .24 .22 .22 .22;
   .27 .21 .20 .18 .17 .17 .17 .17;
   .16 .15 .14 .12 .11 .11 .11 .11;
   .12 .11 .10 .08 .08 .08 .08 .08;
   .09 .08 .07 .06 .06 .06 .06 .05;
    .06 .05 .03 .02 .02 .02 .02 .02];

ns = zeros(NumSamplesLevel,240);
LN = zeros(NumSamplesLevel,1);
N_mat = zeros(NumSamplesLevel,1);
Spec_N = zeros(1,240);

for l = 1:NumSamplesLevel
    N = 0;
    z1 = 0; % critical band rate starts at 0
    n1 = 0; % loudness level starts at 0
    j = 18;
    iz = 1;
    z = 0.1;
    
    for i = 1:21

    % Determines where to start on the slope
       ig = i-1;
       if ig >7;
           ig=7;
       end
       control=1;
       while (z1 < ZUP(i)) || (control==1) % ZUP is the upper limit of the approximated critical band

    % Determines which of the slopes to use
          if n1 < CoreL(l,i),      % Nm is the main loudness level
             j=1;
             while RNS(j) > CoreL(l,i), % the value of j is used below to build a slope
                j=j+1; % j becomes the index at which Nm(i) is first greater than RNS
             end 
          end

    % The flat portions of the loudness graph
         if n1 <= CoreL(l,i),
             z2 = ZUP(i); % z2 becomes the upper limit of the critical band
             n2 = CoreL(l,i);
             N = N + n2*(z2-z1); % Sums the output (N_entire)
             for k = z:0.1:z2      % k goes from z to upper limit of the critical band in steps of 0.1
                ns(l,iz) = n2; % ns is the output, and equals the value of Nm
                if k < (z2-0.05), 
                   iz = iz + 1;
                end
             end
             z = k; % z becomes the last value of k
             z = round(z*10)*0.1; 
         end

    % The sloped portions of the loudness graph
          if n1 > CoreL(l,i),
              n2 = RNS(j);
              if n2 < CoreL(l,i);
                  n2 = CoreL(l,i);
              end
              dz = (n1-n2)/USL(j,ig); % USL = slopes
              dz = round(dz*10)*0.1;
              if dz == 0;
                  dz = 0.1;
              end
              z2 = z1 + dz;
              if z2 > ZUP(i),
                 z2 = ZUP(i);
                 dz = z2-z1;
                 n2 = n1 - dz*USL(j,ig); %USL = slopes
              end
              N = N + dz*(n1+n2)/2; % Sums the output (N_entire)
              for k = z:0.1:z2
                ns(l,iz) = n1 - (k-z1)*USL(j,ig); % ns is the output, USL = slopes
                if k < (z2-0.05),
                   iz = iz + 1;
                end
              end
              z = k;
              z = round(z*10)*0.1;
          end
           if n2 == RNS(j);
               j=j+1;
           end
           if j > 18;
               j = 18;
           end
       n1 = n2;
       z1 = z2;
       z1 = round(z1*10)*0.1;
       control = control+1;
       end
    end

    if N < 0;
        N = 0;
    end

    if N <= 16;
        N = (N*1000+.5)/1000;
    else
        N = (N*100+.5)/100;
    end

    LN(l) = 40*(N + .0005)^.35;

    if LN(l) < 3;
        LN(l) = 3;
    end

    if N >= 1;
        LN(l) = 10*log10(N)/log10(2) + 40;
    end

    N_mat(l) = N;
end

% specific Loudness as a function of Bark number 
for i = 1:240
    Spec_N(i) = max(ns(:,i));
end


% **********************************************************************
% STEP 10 - Apply Temporal Weighting to Arbitrary signals
% **********************************************************************

if method == 1
    
    Loudness_t1 = zeros(NumSamplesLevel,1);
    Loudness_t2 = zeros(NumSamplesLevel,1);
    Loudness = zeros(NumSamplesLevel,1);
    % 1st order low-pass A
    Tau = 3.5e-3;
    A1 = exp(-1 / (SampleRateLevel * DecFactorLevel * Tau));
    B0 = 1 - A1;
    Y1 = 0;
    for i = 1:NumSamplesLevel
        X0 = N_mat(i);
        Y1 = B0 * X0 + A1 * Y1;
        Loudness_t1(i) = Y1;
        if i < NumSamplesLevel - 1
            Xd = (N_mat(i) - X0) / DecFactorLevel;
            for j = 1:DecFactorLevel
                X0 = X0 + Xd;
                Y1 = B0 * X0 + A1 * Y1;
            end
        end
    end
    
    % 1st order low-pass B
    Tau = 70e-3;
    A1 = exp(-1 / (SampleRateLevel * DecFactorLevel * Tau));
    B0 = 1 - A1;
    Y1 = 0;
    for i = 1:NumSamplesLevel
        X0 = N_mat(i);
        Y1 = B0 * X0 + A1 * Y1;
        Loudness_t2(i) = Y1;
        if i < NumSamplesLevel - 1
            Xd = (N_mat(i) - X0) / DecFactorLevel;
            for j = 1:DecFactorLevel
                X0 = X0 + Xd;
                Y1 = B0 * X0 + A1 * Y1;
            end
        end
    end
    
    % combine the filters
    
    for i = 1:NumSamplesLevel
            Loudness(i) = (0.47 * Loudness_t1(i)) + (0.53 * Loudness_t2(i));
    end

    % Decimate signal for decreased computation time by factor of fs =
    % 2 Hz 
    
    Total_Loudness = decimate(Loudness,DecFactorLoudness);
    
    ns_dec = zeros(NumSamplesLoudness,240);
    for i = 1:240
        ns_dec(:,i) = decimate(ns(:,i),DecFactorLoudness);
    end

    % statistics
    Nmax = max(Total_Loudness); % result presented by the standard
    N5 = prctile(Total_Loudness,95); % result presented by the standard
    Nmean = mean(Total_Loudness);
    Nstd = std(Total_Loudness);
    N1 = prctile(Total_Loudness,99);
    N2 = prctile(Total_Loudness,98);
    N3 = prctile(Total_Loudness,97);
    N4 = prctile(Total_Loudness,96);
    N10 = prctile(Total_Loudness,90);
    N20 = prctile(Total_Loudness,80);
    N30 = prctile(Total_Loudness,70);
    N40 = prctile(Total_Loudness,60);
    N50 = median(Total_Loudness);
    N60 = prctile(Total_Loudness,40);
    N70 = prctile(Total_Loudness,30);
    N80 = prctile(Total_Loudness,20);
    N90 = prctile(Total_Loudness,10);
    Nmin = min(Total_Loudness);
    
    data_Total_Loudness = [Nmean;Nstd;Nmax;N1;N2;N3;N4;N5;N10;N20;N30;N40;N50;N60;N70;N80;N90;Nmin];

    
    tPeriod = 2e-3; % 2 ms
    timePoints = (0:length(Total_Loudness)-1)' * tPeriod;
    
    fig1 = figure('Name','Time-varying Dynamic Loudness (C&F) and Sharpness Statistics');
    table1 = uitable('Data',data_Total_Loudness,...
        'ColumnName',{'Loudness'},...
        'RowName',{'Mean','Standard deviation','Maximum',...
        'N1','N2','N3','N4',...
        'N5','N10','N20','N30','N40','N50 (median)','N60',...
        'N70','N80','N90','Minimum'});
    
    [~,tables] = disptables(fig1,table1); % AARAE function
    
    OUT.tables = tables;
    
    % Figure for charts
    figure('Name',['Time-varying method Loudness (ISO532-1) of ',name])
    set(gca,'FontSize',14.);
    subplot(2,2,1:2)
    % Time-varying loudness 
    plot(timePoints,Total_Loudness,'r-');
    title ('Time-Varying Loudness');
    xlabel('Time (s)');
    ylabel('Loudness (sone)','Color','k');
    hold off;
    
    subplot(2,2,4)
    % Time-averaged specific loudness as a fucntion of critical band
    plot((1:240)/10,Spec_N,'r-');
    ax=gca;
    ax.Title.String = 'Time-Averaged Specific Loudness';
    ax.XLabel.String = 'Critical Band Rate (Bark)';
    ax.XLim = [0 length(Spec_N)/10+1];
    % ax.XTickLabel = {'0','5','10','15','20','25'};
    ax.YLabel.String = 'Loudness (sones/Bark)';
    hold off;

    subplot(2,2,3)
     % Specific loudness spectrogram
    imagesc(timePoints, 0.1:0.1:24, ns_dec');
    % cH = colorbar;
    set(gca,'YDir','normal');
    ax=gca;
    axis tight;
    ax.Title.String = 'Specific Loudness';
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Critical Band Rate (Bark)';
    hold off;
    
elseif method == 0
    data = [N;LN];
    fig1 = figure('Name','Overall Loudness and Loudness Level');
    table1 = uitable('Data',data,...
        'RowName',{'N (sone)','LN (phone)'});
    
    [~,tables] = disptables(fig1,table1); % AARAE function
    
    OUT.tables = tables;
  
    % Figure for charts
    figure('Name',['Stationary method Loudness (ISO532-1) of ',name])
    % Time-averaged specific loudness as a fucntion of critical band
    plot((1:240)/10,Spec_N,'r-');
    ax=gca;
    ax.FontSize = 14.;
    ax.Title.String = 'Time-Averaged Specific Loudness';
    ax.XLabel.String = 'Critical Band Rate (Bark)';
    ax.XLim = [0 length(Spec_N)/10+1];
    ax.YLabel.String = 'Loudness (sones/Bark)';
    hold off;
    
end

% ******** AARAE RESULTS LEAVES **********

if isstruct(IN)
    if method == 0
        % Time-averaged Specific loudness results leaf
        doresultleaf(Spec_N','Specific Loudness [sones/Bark]',{'Critical Band Rate'},...
            'Critical Band Rate',[0.1:0.1:24]','Bark',true,...
            'loudnesstype', {'Loudness over critical band'}, 'categorical', [],...
            'name','Time_averaged_specific_loudness');
    
    elseif method == 1  
        % Time-varying loudness results leaf
        doresultleaf(Total_Loudness,'Loudness [sone]',{'time'},...
            'Time',timePoints','s',true,...
            'loudnesstype', {'Loudness over time'}, 'categorical',[],...
            'name','Time_varying_loudness');
        
        % Time-averaged Specific loudness results leaf
        doresultleaf(Spec_N','Specific Loudness [sones/Bark]',{'Critical Band Rate'},...
            'Critical Band Rate',[0.1:0.1:24]','Bark',true,...
            'loudnesstype', {'Loudness over critical band'}, 'categorical', [],...
            'name','Time_averaged_specific_loudness');
        
        % Specific loudness spectrogram results leaf
        doresultleaf(ns_dec','Specific Loudness [sones/Bark]',{'time','Critical Band Rate'},...
            'Critical Band Rate',[0.1:0.1:24]','Bark',true,...
            'Time',timePoints','s',true,...
            'loudnesstype', {'Loudness over critical band'}, 'categorical', [],...
            'name','Time_varying_specific_loudness');
    end
end

OUT.funcallback.name = 'Loudness_ISO532_1.m';
OUT.funcallback.inarg = {field, method, TIME_SKIP};

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