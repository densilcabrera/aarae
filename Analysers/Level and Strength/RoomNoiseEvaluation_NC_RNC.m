function out = RoomNoiseEvaluation_NC_RNC(in, fs, cal, showpercentiles)
% This function calculates background noise ratings based on ANSI S12.2
% (2008), namely Noise Criterion (NC) and Room Noise Criterion (RNC).
%
% code by Guy Hopkins, Nicholas Lynar & Densil Cabrera
% version 0 (not validated yet)

if isstruct(in)
    in = choose_from_higher_dimensions(in,1,1); 
    audio = in.audio;
    fs = in.fs;
    if isfield(in,'cal')
        cal = in.cal;
    else
        cal = 0;
        disp('This audio signal has not been calibrated.')
    end
else
    audio = in;
    if nargin < 3
        cal = inputdlg({'Calibration offset (dB)'},...
            'cal',1,{'0'});
        cal = str2num(char(cal));
    end
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
            'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

if nargin < 4,
    showpercentiles = 1;
    param = inputdlg({'Calibration offset (dB)'; ...
        'Show percentiles [0 | 1]'}, ...
        'Analysis and display parameters',1, ...
        {num2str(cal); ...
        num2str(showpercentiles)});
    param = str2num(char(param));
    if length(param) < 2, param = []; end
    if ~isempty(param)
        
        cal = param(1);
        showpercentiles = param(2);
    end
end

if ~isempty(audio) && ~isempty(fs) && ~isempty(cal)
    [len,chans,bands] = size(audio);
    
    % mix down multiband input
    if bands > 1
        audio = sum(audio,3);
        disp('Multiband audio has been mixed')
    end
    
    % calibration adjustment, and mix down to 1 channel
    if chans > 1
        if length(cal) == chans
            audio = audio .* 10.^(repmat(cal,[len,1])./20);
        elseif length(cal) == 1
            audio = audio .* 10.^(cal./20);
        else
            % give up!
            disp('Calibration vector does not clearly map to channels')
            out = [];
            return
        end
        audio = mean(audio,2);
        disp('Multiple audio channels have been averaged')
    else
        % calibration adjustment for single channel input
        audio = audio .* 10.^(cal./20);
    end
    % Octave band filterbank from AARAE: in Processors/Filterbanks
    frequencies = [16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000];
    audiooct = octbandfilter_viaFFT(audio,fs,frequencies);
    
    % A-weight - use AARAE's A-weighting filter: in Processors/Filters
    audioA = Aweight(audio,fs);
    
    % CALCULATE TEMPORAL INTEGRATION
    % FILTER DESIGN
    tau = 0.125; % fast time constant
    E = exp(-1/(tau*fs)); % exponential term
    b = 1 - E; % filter numerator (adjusts gain to compensate for denominator)
    a = [1 -E];% filter denominator
    % rectify, integrate and convert to decibels
    
    Ifast=filter(b,a,abs(audiooct)).^2;
    IfastA=filter(b,a,abs(audioA)).^2;
    
    tau = 1; % slow time constant
    E = exp(-1/(tau*fs)); % exponential term
    b = 1 - E; % filter numerator (adjusts gain to compensate for denominator)
    a = [1 -E];% filter denominator
    % rectify, integrate and convert to decibels
    
    Islow=filter(b,a,abs(audiooct)).^2;
    IslowA=filter(b,a,abs(audioA)).^2;
    
    %SIL Calculations first method outlined for NC rating
    
    L27= 10*log10(mean(Ifast(:,1,6))); %500Hz
    L30= 10*log10(mean(Ifast(:,1,7))); %1000Hz
    L33= 10*log10(mean(Ifast(:,1,8))); %2000Hz
    L36= 10*log10(mean(Ifast(:,1,9))); %4000Hz
    SIL= 0.25 .*(L27+L30+L33+L36); %Average of the Octave bands 500-4000Hz
    out.SIL = round(SIL);%Sound Interference Level
    %Noise Criterion curves
    NCcurves = ...
   [90.0	90.0	84.0	79.0	75.0	72.0	71.0	70.0	68.0	68.0
90.0	89.6	83.2	78.2	74.2	71.2	69.8	68.8	67.0	66.8
90.0	89.2	82.4	77.4	73.4	70.4	68.6	67.6	66.0	65.6
90.0	88.8	81.6	76.6	72.6	69.6	67.4	66.4	65.0	64.4
90.0	88.4	80.8	75.8	71.8	68.8	66.2	65.2	64.0	63.2
90.0	88.0	80.0	75.0	71.0	68.0	65.0	64.0	63.0	62.0
90.0	87.4	79.4	74.2	70.0	67.0	64.0	63.0	62.0	61.0
90.0	86.8	78.8	73.4	69.0	66.0	63.0	62.0	61.0	60.0
90.0	86.2	78.2	72.6	68.0	65.0	62.0	61.0	60.0	59.0
90.0	85.6	77.6	71.8	67.0	64.0	61.0	60.0	59.0	58.0
90.0	85.0	77.0	71.0	66.0	63.0	60.0	59.0	58.0	57.0
89.8	84.4	76.4	70.2	65.2	62.0	59.2	58.0	57.0	56.0
89.6	83.8	75.8	69.4	64.4	61.0	58.4	57.0	56.0	55.0
89.4	83.2	75.2	68.6	63.6	60.0	57.6	56.0	55.0	54.0
89.2	82.6	74.6	67.8	62.8	59.0	56.8	55.0	54.0	53.0
89.0	82.0	74.0	67.0	62.0	58.0	56.0	54.0	53.0	52.0
88.6	81.4	73.4	66.4	61.2	57.2	55.0	53.0	52.0	51.0
88.2	80.8	72.8	65.8	60.4	56.4	54.0	52.0	51.0	50.0
87.8	80.2	72.2	65.2	59.6	55.6	53.0	51.0	50.0	49.0
87.4	79.6	71.6	64.6	58.8	54.8	52.0	50.0	49.0	48.0
87.0	79.0	71.0	64.0	58.0	54.0	51.0	49.0	48.0	47.0
86.6	78.4	70.2	63.2	57.2	53.0	50.0	48.0	47.0	46.0
86.2	77.8	69.4	62.4	56.4	52.0	49.0	47.0	46.0	45.0
85.8	77.2	68.6	61.6	55.6	51.0	48.0	46.0	45.0	44.0
85.4	76.6	67.8	60.8	54.8	50.0	47.0	45.0	44.0	43.0
85.0	76.0	67.0	60.0	54.0	49.0	46.0	44.0	43.0	42.0
84.8	75.6	66.4	59.2	53.2	48.0	45.0	43.0	42.0	41.0
84.6	75.2	65.8	58.4	52.4	47.0	44.0	42.0	41.0	40.0
84.4	74.8	65.2	57.6	51.6	46.0	43.0	41.0	40.0	39.0
84.2	74.4	64.6	56.8	50.8	45.0	42.0	40.0	39.0	38.0
84.0	74.0	64.0	56.0	50.0	44.0	41.0	39.0	38.0	37.0
83.6	73.4	63.2	55.2	49.0	43.2	40.0	38.0	37.0	36.0
83.2	72.8	62.4	54.4	48.0	42.4	39.0	37.0	36.0	35.0
82.8	72.2	61.6	53.6	47.0	41.6	38.0	36.0	35.0	34.0
82.4	71.6	60.8	52.8	46.0	40.8	37.0	35.0	34.0	33.0
82.0	71.0	60.0	52.0	45.0	40.0	36.0	34.0	33.0	32.0
81.8	70.4	59.4	51.2	44.2	39.0	35.2	33.0	32.0	31.0
81.6	69.8	58.8	50.4	43.4	38.0	34.4	32.0	31.0	30.0
81.4	69.2	58.2	49.6	42.6	37.0	33.6	31.0	30.0	29.0
81.2	68.6	57.6	48.8	41.8	36.0	32.8	30.0	29.0	28.0
81.0	68.0	57.0	48.0	41.0	35.0	32.0	29.0	28.0	27.0
80.8	67.4	56.4	47.2	40.2	34.2	31.0	28.0	26.8	26.0
80.6	66.8	55.8	46.4	39.4	33.4	30.0	27.0	25.6	25.0
80.4	66.2	55.2	45.6	38.6	32.6	29.0	26.0	24.4	24.0
80.2	65.6	54.6	44.8	37.8	31.8	28.0	25.0	23.2	23.0
80.0	65.0	54.0	44.0	37.0	31.0	27.0	24.0	22.0	22.0
79.8	64.6	53.2	43.2	36.2	30.0	26.0	23.2	21.0	20.8
79.6	64.2	52.4	42.4	35.4	29.0	25.0	22.4	20.0	19.6
79.4	63.8	51.6	41.6	34.6	28.0	24.0	21.6	19.0	18.4
79.2	63.4	50.8	40.8	33.8	27.0	23.0	20.8	18.0	17.2
79.0	63.0	50.0	40.0	33.0	26.0	22.0	20.0	17.0	16.0
78.8	62.6	49.4	39.2	32.0	25.2	21.2	18.8	16.0	15.0
78.6	62.2	48.8	38.4	31.0	24.4	20.4	17.6	15.0	14.0
78.4	61.8	48.2	37.6	30.0	23.6	19.6	16.4	14.0	13.0
78.2	61.4	47.6	36.8	29.0	22.8	18.8	15.2	13.0	12.0
78.0	61.0	47.0	36.0	28.0	22.0	18.0	14.0	12.0	11.0];

spectrum = 10*log10(mean(Ifast));
spectrum = permute(spectrum,[1,3,2]);
%Finds appropiate curve using the tangency method
NCoctind=zeros(10,1);
for k = 1:10
    ind= find(NCcurves(:,k)>=spectrum(k),1,'last');
    if ~isempty(ind)
        NCoctind(k) = ind;
    end
end
if min(NCoctind)~=0
    NCval = 70:-1:15;
    NCoct=NCval(:,NCoctind);
    maxNCoctind=min(NCoctind);

    out.NC = max(NCoct);
else
    maxNCoctind = [];
    out.NC = [];
    disp('NC out of range')
end

%NCSILind = (70-(out.SIL-15));
NCSILind=70-out.SIL-1;
if ~(NCSILind<1) && ~(NCSILind>56)
    NCSILcurve= NCcurves(NCSILind,:);
else
    NCSILcurve=[];
    disp('NC out of range')
end



%if the measured spectrum in any octave band does not exceed any of the octave bands of that NC(SIL) the SIL method may be used. Otherwise Tangency.    
    % Weighted sum of lowest three bands
    RNClowbands = ((Ifast(:,1,1))-((10^(14/20))*(2*10^-5))) ... % 16 Hz
        + Ifast(:,1,2) ... % 31.5 Hz
        + ((Ifast(:,1,3))+((10^(14/20))*(2*10^-5))); % 63 Hz
    % 125 Hz octave band
    RNC125 = Ifast (:,1,4);
    % concatenate with 125 Hz octave band
    % RNClowbands = cat(3,RNClowbands,Ifast(:,1,4));
   
    
    % calculate SPL of 100 ms windows
    nwindows = floor(len ./(fs*0.1));
    winlen = floor(fs*0.1);
    [RNClowoctIntensity,RNC125LevelIntensity] = deal(zeros(nwindows,1));
    for n = 1:nwindows
        RNClowoctIntensity =  ...
            RNClowbands((1):winlen:((n-1)*winlen+winlen));
        RNC125LevelIntensity =  ...
            RNC125((1):winlen:((n-1)*winlen+winlen));
            %RNC125(((n-1)*winlen+1):((n-1)*winlen+winlen));
        
    end
            t = 0.1*((1:nwindows)-1);
            doresultleaf(10*log10([RNClowoctIntensity,RNC125LevelIntensity]), ['LowFreqBandLevels' '[dB]'],{'time'},...
                 'time',     t,      's',           true,...
                 'bands', {'low' '125 Hz'}, 'categorical', [],...
                 'name','RNC_LF_Fluctuation');
    
    
    %Lowband Levels
    lowRNCLeq = 10*log10(mean(RNClowoctIntensity));%time-averaged level of lowband spectra
    lowRNCLmean = mean(10*log10(RNClowoctIntensity)); %Mean Lp
    lowRNCL10 = 10*log10(prctile(RNClowoctIntensity,90));%L10
    lowRNCLmax = 10*log10(max(RNClowoctIntensity));%Lmax
    lowRNCLStDev= std(10*log10(RNClowoctIntensity));%standard deviation
    
    %125 Hz band levels
    
    RNC125Leq = 10*log10(mean(RNC125LevelIntensity)); %time-averaged level of 125Hz
    RNC125Lmean = mean(10*log10(RNC125LevelIntensity)); %mean Lp
    RNC125L10 = 10*log10(prctile(RNC125LevelIntensity,90));%L10
    RNC125Lmax = 10*log10(max(RNC125LevelIntensity));%Lmax
    RNC125LStDev= std(10*log10(RNC125LevelIntensity));%standard deviation
    
    %determination of large random fluctuations
    LowMaxMinusLeq=lowRNCLmax-lowRNCLeq;
    LowL10minusLeq=lowRNCL10-lowRNCLeq;
    MaxMinus125Leq=RNC125Lmax-RNC125Leq;
    L10minus125Leq=RNC125L10-RNC125Leq;
    
    %correction factor 31.5 Hz (lowbands) large random fluctuations present
        N = nwindows;            %number of 100ms windows in measured signal.
        Lajlow=2.*((10*log10(RNClowoctIntensity))-lowRNCLmean);
        Klfalow= 10.*log10((1/N).*(sum(10.^(Lajlow/10))));
        corr_factorlow = Klfalow-(lowRNCLeq-lowRNCLmean); %correction factor for 3 band sum, centered at 31.5 Hz.
        %if surging present apply this above correction factor!!!
        
        %correction factor 31.5 Hz (lowbands) no strong surging present
        
        corr_factor_lowsurge=0.46.*((lowRNCLStDev)^2);
        
         %correction factor 125 Hz (lowbands) Surging present
        N = nwindows;            %%number of 100ms windows in measured signal.
        Laj125= 1.25 .*((10*log10(RNC125LevelIntensity))-RNC125Lmean);
        Klfa125= 10.*log10((1/N).*(sum(10.^(Laj125/10))));
        corr_factor125 = Klfa125-(RNC125Leq-RNC125Lmean); %correction factor for 3 band sum, placed in 31.5 hz band.
        %if surging present apply this above correction factor!!!
        %correction factor 125 Hz (lowbands) no strong surging present
        
        corr_factor_l25surge=0.18.*((RNC125LStDev)^2);
        
    %Determination of correction factor. Derived from 'Screening criteria
    %to determine presence of large random fluctuations.'
    
    if LowMaxMinusLeq >= 7;
        out.corrfactorlow= corr_factorlow;
    elseif LowMaxMinusLeq < 7
        out.corrfactorlow = corr_factor_lowsurge;
    elseif LowL10minusLeq >=3.5;
       out.corrfactorlow=corr_factorlow;
    elseif LowL10minusLeq <3.5;
       out.corrfactorlow=corr_factor_lowsurge;
    end
    if MaxMinus125Leq >=6;
       out.corrfactor125=corr_factor125;
    elseif MaxMinus125Leq <6;
       out.corrfactor125=corr_factor_l25surge;
    elseif L10minus125Leq >= 3;
        out.corrfactor125=corr_factor125;
    elseif L10minus125Leq < 3;
        out.corrfactor125=corr_factor_l25surge;
    end
    
    % A-weighted Leq
    out.LAeqslow = 10*log10(mean(IslowA));
    out.LAeqfast = 10*log10(mean(IfastA));
    
    % Octave band Leq (unweighted), and other octave band values below
    %outputs required in table: mean level, Standard deviation, LEQ,
    %Leq-mean level, Corrections to leq, final spectrum,
    
    %single values: LAEq, NC rating, RNC rating. (method determination)
    out.Leq = 10*log10(mean(Ifast)); %octave band Leq, from measurement (fast)
    out.Leq = permute(out.Leq,[1,3,2]); %(Ifast)
    out.StandardDev = std(10*log10(Ifast));%standard deviation at each octave band
    out.StandardDev = permute(out.StandardDev, [1,3,2]);
    out.Corrections = [0 out.corrfactorlow 0 out.corrfactor125 0 0 0 0 0 0]; %correction factors at 31.5Hz and 125 hz
    out.Finalspectrum = out.Leq+out.Corrections;%final spectrum with correction added
    out.Finalspectrum = out.Finalspectrum .* [-inf 1 -inf 1 1 1 1 1 1 1];
    out.Leqslow = 10*log10(mean(Islow));
    out.Leqfast = 10*log10(mean(Ifast));
    out.Leqslow = permute(out.Leqslow,[1,3,2]);
    out.Leqfast = permute(out.Leqfast,[1,3,2]);
    
    %RNC Curves
    RNCcurves = ...
[101.00	96.00	91.00	86.00	81.00	76.00	72.00	68.00	64.00	60.00
100.00	95.00	90.00	85.00	80.00	75.00	71.00	67.00	63.00	59.00
99.00	94.00	89.00	84.00	79.00	74.00	70.00	66.00	62.00	58.00
98.00	93.00	88.00	83.00	78.00	73.00	69.00	65.00	61.00	57.00
97.00	92.00	87.00	82.00	77.00	72.00	68.00	64.00	60.00	56.00
96.00	91.00	86.00	81.00	76.00	71.00	67.00	63.00	59.00	55.00
95.00	90.00	85.00	80.00	75.00	70.00	66.00	62.00	58.00	54.00
94.00	89.00	84.00	79.00	74.00	69.00	65.00	61.00	57.00	53.00
93.00	88.00	83.00	78.00	73.00	68.00	64.00	60.00	56.00	52.00
92.00	87.00	82.00	77.00	72.00	67.00	63.00	59.00	55.00	51.00
91.00	86.00	81.00	76.00	71.00	66.00	62.00	58.00	54.00	50.00
90.00	85.00	80.00	75.00	70.00	65.00	61.00	57.00	53.00	49.00
89.00	84.00	79.00	74.00	69.00	64.00	60.00	56.00	52.00	48.00
88.00	83.00	78.00	73.00	68.00	63.00	59.00	55.00	51.00	47.00
87.00	82.00	77.00	72.00	67.00	62.00	58.00	54.00	50.00	46.00
86.00	81.00	76.00	71.00	66.00	61.00	57.00	53.00	49.00	45.00
85.00	80.00	75.00	70.00	65.00	60.00	56.00	52.00	48.00	44.00
84.00	79.00	74.00	69.00	64.00	59.00	55.00	51.00	47.00	43.00
83.00	78.00	73.00	68.00	63.00	58.00	54.00	50.00	46.00	42.00
82.00	77.00	72.00	67.00	62.00	57.00	53.00	49.00	45.00	41.00
81.00	76.00	71.00	66.00	61.00	56.00	52.00	48.00	44.00	40.00
80.67	75.50	70.33	65.17	60.00	55.00	51.00	47.00	43.00	39.00
80.33	75.00	69.67	64.33	59.00	54.00	50.00	46.00	42.00	38.00
80.00	74.50	69.00	63.50	58.00	53.00	49.00	45.00	41.00	37.00
79.67	74.00	68.33	62.67	57.00	52.00	48.00	44.00	40.00	36.00
79.33	73.50	67.67	61.83	56.00	51.00	47.00	43.00	39.00	35.00
79.00	73.00	67.00	61.00	55.00	50.00	46.00	42.00	38.00	34.00
78.67	72.50	66.33	60.17	54.00	49.00	45.00	41.00	37.00	33.00
78.33	72.00	65.67	59.33	53.00	48.00	44.00	40.00	36.00	32.00
78.00	71.50	65.00	58.50	52.00	47.00	43.00	39.00	35.00	31.00
77.67	71.00	64.33	57.67	51.00	46.00	42.00	38.00	34.00	30.00
77.33	70.50	63.67	56.83	50.00	45.00	41.00	37.00	33.00	29.00
77.00	70.00	63.00	56.00	49.00	44.00	40.00	36.00	32.00	28.00
76.67	69.50	62.33	55.17	48.00	43.00	39.00	35.00	31.00	27.00
76.33	69.00	61.67	54.33	47.00	42.00	38.00	34.00	30.00	26.00
76.00	68.50	61.00	53.50	46.00	41.00	37.00	33.00	29.00	25.00
75.67	68.00	60.33	52.67	45.00	40.00	36.00	32.00	28.00	24.00
75.33	67.50	59.67	51.83	44.00	39.00	35.00	31.00	27.00	23.00
75.00	67.00	59.00	51.00	43.00	38.00	34.00	30.00	26.00	22.00
74.67	66.50	58.33	50.17	42.00	37.00	33.00	29.00	25.00	21.00
74.33	66.00	57.67	49.33	41.00	36.00	32.00	28.00	24.00	20.00
74.00	65.50	57.00	48.50	40.00	35.00	31.00	27.00	23.00	19.00
73.67	65.00	56.33	47.67	39.00	34.00	30.00	26.00	22.00	18.00
73.33	64.50	55.67	46.83	38.00	33.00	29.00	25.00	21.00	17.00
73.00	64.00	55.00	46.00	37.00	32.00	28.00	24.00	20.00	16.00
72.67	63.50	54.33	45.17	36.00	31.00	27.00	23.00	19.00	15.00
72.33	63.00	53.67	44.33	35.00	30.00	26.00	22.00	18.00	14.00
72.00	62.50	53.00	43.50	34.00	29.00	25.00	21.00	17.00	13.00
71.67	62.00	52.33	42.67	33.00	28.00	24.00	20.00	16.00	12.00
71.33	61.50	51.67	41.83	32.00	27.00	23.00	19.00	15.00	11.00
71.00	61.00	51.00	41.00	31.00	26.00	22.00	18.00	14.00	10.00
70.67	60.50	50.33	40.17	30.00	25.00	21.00	17.00	13.00	9.00
70.33	60.00	49.67	39.33	29.00	24.00	20.00	16.00	12.00	8.00
70.00	59.50	49.00	38.50	28.00	23.00	19.00	15.00	11.00	7.00
69.67	59.00	48.33	37.67	27.00	22.00	18.00	14.00	10.00	6.00
69.33	58.50	47.67	36.83	26.00	21.00	17.00	13.00	9.00	5.00
69.00	58.00	47.00	36.00	25.00	20.00	16.00	12.00	8.00	4.00
68.67	57.50	46.33	35.17	24.00	19.00	15.00	11.00	7.00	3.00
68.33	57.00	45.67	34.33	23.00	18.00	14.00	10.00	6.00	2.00
68.00	56.50	45.00	33.50	22.00	17.00	13.00	9.00	5.00	1.00
67.67	56.00	44.33	32.67	21.00	16.00	12.00	8.00	4.00	0.00];
%Finds appropiate curve using the tangency method
spectrumrnc = out.Finalspectrum;
RNCoctind=zeros(10,1);
for k = 1:10
    ind= find(RNCcurves(:,k)>=spectrumrnc(k),1,'last');
    if ~isempty(ind)
        RNCoctind(k) = ind;
    end
end
if min(RNCoctind)~=0
    RNCval = 70:-1:10;
    RNCoct=RNCval(:,RNCoctind);
    maxRNCoctind=min(RNCoctind);
    out.RNC = max(RNCoct);
else
    maxRNCoctind = [];
    out.RNC = [];
    disp('RNC out of range')
end


    
    out.Lmaxslow = 10*log10(max(Islow));
    out.Lmaxfast = 10*log10(max(Ifast));
    out.Lmaxslow = permute(out.Lmaxslow,[1,3,2]);
    out.Lmaxfast = permute(out.Lmaxfast,[1,3,2]);
    
    out.L10slow = 10*log10(prctile(Islow,90));
    out.L10fast = 10*log10(prctile(Ifast,90));
    out.L10slow = permute(out.L10slow,[1,3,2]);
    out.L10fast = permute(out.L10fast,[1,3,2]);
    
    out.L90slow = 10*log10(prctile(Islow,10));
    out.L90fast = 10*log10(prctile(Ifast,10));
    out.L90slow = permute(out.L90slow,[1,3,2]);
    out.L90fast = permute(out.L90fast,[1,3,2]);
    
    out.funcallback.name = 'RoomNoiseEvalueation_NC_RNC.m';
    out.funcallback.inarg = {fs, cal, showpercentiles};
    
    
    
    % OUTPUT TABLE
    fig1 = figure('Name','ANSI 12.2 (2008) values');
    table1 = uitable('Data',...
        [out.Leq;out.StandardDev;out.Corrections;...
        out.Finalspectrum],...
        'ColumnName',num2cell(frequencies),...    
        'RowName',{'Leq (fast)',...
        'Standard Deviation','Correction Factor','Final Spectrum'});
    
    
    % Octave band values table - fast
    table3 = uitable('Data',...
      [LowMaxMinusLeq MaxMinus125Leq;LowL10minusLeq L10minus125Leq],...
        'ColumnName',num2cell([31.5 125]),...    
        'RowName',{'Max minus Leq',...
        'L10 minus Leq',});


     % Single number values table
    table2 = uitable('Data',...
        [out.LAeqslow;...
        out.SIL;out.NC;out.RNC],...
        'ColumnName', {'Value'},...
        'RowName',{'LAeq (slow)';...
        'NC(SIL)';'NC(tangency)';'RNC(tangency)'});
    [~,tables] = disptables(fig1,[table1 table2 table3]);
    out.tables = tables;
    
    
    
    
    figure('name', 'Octave Band Spectrum: NC and RNC')
    
    %['Octave Band Spectrum, fast, LAeq = ',num2str(out.LAeqfast),' dBA']
    width = 0.5;
    subplot(2,1,1)
    ymax = 10*ceil(max(max(out.Lmaxfast+5))/10);
    ymin = 10*floor(min(min(out.L90fast))/10);
    title(['NC ',num2str(out.NC)])
    bar(1:length(frequencies),out.Leqfast',width,'FaceColor',[1,0.3,0.3],...
        'EdgeColor',[0,0,0],'DisplayName', 'Leq','BaseValue',ymin);
    hold on
    title(['NC ',num2str(out.NC)])
    %title(['Fast LAeq = ',num2str(out.LAeqfast),' dBA'])
    % x-axis
    set(gca,'XTickLabel',num2cell(frequencies))
    if (frequencies(1) ~= 1) && (frequencies(end) ~= length(frequencies))
        xlabel('Octave Band Centre Frequency (Hz)')
    else
        xlabel('Band')
    end
    
    % y-axis
    ylabel('Level (dB)')
    ylim([ymin ymax])
    if ~isempty(maxNCoctind)
    plot(1:length(frequencies),NCcurves(maxNCoctind,:),'Color',[0,0,1],...
        'DisplayName','NC curve')
    else
        text(1,ymax-0.13*(ymax-ymin),'NC out of range')
    end
    if ~isempty(NCSILcurve)
    plot(1:length(frequencies),NCSILcurve,'Color',[0,1,0],...
        'DisplayName','NC(SIL) curve')
    else
         text(1,ymax-0.13*(ymax-ymin),'NC out of range')
    end

    if showpercentiles
    
        
        legend('show','Location','EastOutside');
        hold off
    else
        legend 'off'
    end
    
    for k = 1:length(frequencies)
        text(k-0.25,ymax-(ymax-ymin)*0.04, ...
            num2str(round(out.Leqfast(k)*10)/10),'Color',[1,0.3,0.3])
    end
    
    
    % ***
    subplot(2,1,2)
    bar(1:length(frequencies),out.Finalspectrum',width,'FaceColor',[0.3,0.7,0.3],...
        'EdgeColor',[0,0,0],'DisplayName', 'Leq','BaseValue',ymin);
    hold on
    %title(['Slow LAeq = ',num2str(out.LAeqslow),' dBA'])
    title(['RNC ',num2str(out.RNC)])
    % x-axis
    set(gca,'XTickLabel',num2cell(frequencies))
    if (frequencies(1) ~= 1) && (frequencies(end) ~= length(frequencies))
        xlabel('Octave Band Centre Frequency (Hz)')
    else
        xlabel('Band')
    end
    
    % y-axis
    ylabel('Level (dB)')
    ylim([ymin ymax])
    if ~isempty(maxRNCoctind)
    plot(1:length(frequencies),RNCcurves(maxRNCoctind,:),'Color',[1,0,0],...
        'DisplayName','RNC curve')
    else
        text(1,ymax-0.13*(ymax-ymin),'RNC out of range')
    end
    if LowMaxMinusLeq >= 7;
        text(8,ymax-0.13*(ymax-ymin),'Surge/Fluct.')
    elseif LowL10minusLeq >=3.5;
        text(8,ymax-0.13*(ymax-ymin),'Surge/Fluct.')
    elseif MaxMinus125Leq >=6;
        text(8,ymax-0.13*(ymax-ymin),'Surge/Fluct.')
    elseif L10minus125Leq >= 3;
        text(8,ymax-0.13*(ymax-ymin),'Surge/Fluct.')
    end
    if showpercentiles
        
        
        legend('show','Location','EastOutside');
        hold off
    else
        legend 'off'
    end
    
    for k = 1:length(frequencies)
        text(k-0.25,ymax-(ymax-ymin)*0.04, ...
            num2str(round(out.Finalspectrum(k)*10)/10),'Color',[0.3,0.7,0.3])
    end
    
else
    out = [];
end

%**************************************************************************
% Copyright (c) 2014, Guy Hopkins, Nicholas Lynar & Densil Cabrera
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
%  * Neither the name of the University of Sydney nor the names of its contributors
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