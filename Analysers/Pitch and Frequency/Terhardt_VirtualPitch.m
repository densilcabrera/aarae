function [OUT] = Terhardt_VirtualPitch(IN,timestep,PitchShift,MINWEIGHT,MMAX,SPWEIGHT,maxfrequency,ATune,chromaorder,cal,fs)
% This function analyses the input audio waveform in terms of pitch. It can
% predict multiple simultaneous pitch sensations, including spectral and
% virtual pitches, and their respective pitch strengths (salience).
%
% This function implements the pitch analysis algorithm described in:
%
% E. Terhardt, G. Stoll, and M. Seewann,  1982,
% Algorithm for Extraction of Pitch and Pitch Salience from Complex Tonal
% Signals. Journal of the Acoustical Society of America, 71(3), 679-688.
%
% It also includes some extensions described in:
%
% R. Parncutt. 1989, Harmony: A Psychoacoustical Approach,
% Berlin: Springer-Verlag.
%
% The main code is a port of the implementation in PsySound2 (by Densil
% Cabrera) which was written in Pascal. It was ported to Matlab (and AARAE)
% by Densil Cabrera. (Note that previously this was ported to the PsySound3
% Matlab project, by Matt Flax, but that implementation is no longer
% working.)
%
% Analysis is done in a succession of 80 ms windows (using a Hann window).
% Therefore the time step between windows should be less than 80 ms - this
% can be adjusted by the user.
%
% Version 1.00 May 2015

% ***********************************************************************
if isstruct(IN)
    if isfield(IN,'cal')
        cal = IN.cal;
    else
        cal = 70; % a useable cal offset for typical recordings
    end
end

if nargin ==1

    param = inputdlg({'Time step between windows (ms) - normally < 80 ms';... % These are the input box titles in the
        'Apply pitch shift [0 | 1]';...
        'Minimum weight to output';...
        'Maximum number of subharmonics in pitch calculation';...
        'Spectral weight multiplier';...
        'Maximum frequency for pitch calculation';...
        'Tuning A in Hz';...
        'Calibration offset (dB)';...
        'Chroma in ascending order [0] or cycle-of-fifths order [1]'},...
        'Terhardt Virtual Pitch Analysis Settings',... % This is the dialog window title.
        [1 60],...
        {'20';'1';'0.1';'12';'0.5';'5000';'440';num2str(cal);'1'}); % And the preset answers for your dialog.
    
    param = str2num(char(param));
    
    if length(param) < 1, param = []; end
    if ~isempty(param)
        timestep = param(1);
        PitchShift = param(2);
        if PitchShift ~= (0 | 1), PitchShift = 0; end
        MINWEIGHT = param(3);
        MMAX = param(4);
        SPWEIGHT = param(5);
        maxfrequency = param(6);
        ATune = param(7);
        cal = param(8);
        chromaorder = param(9);
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
    IN = choose_from_higher_dimensions(IN,3,1); % maximum 3 input audio dimensions
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
elseif ~isempty(param) || nargin > 1
    audio = IN;
    % fs and cal are required inputs for standalone
end
% *************************************************************************





% Check that the required data exists for analysis to run
if ~isempty(audio) && ~isempty(fs) && ~isempty(cal)
    
    
    % the audio's length, number of channels, and number of bands
    [len,chans,bands] = size(audio);
    
    if bands > 1
        audio = sum(audio,3); % mixdown bands if multiband
    end
    
    if chans > 1
        audio = mean(audio,2); % mixdown channels
    end
    
    audio = cal_reset_aarae(audio,0,cal);
    
    % Most of the following are set through input arguments or dialog box
    %PitchShift = 1;
    MAXTONES = 128;	%{Maximum number of tones in pitch calculation}
    %MINWEIGHT = 0.1;%{Minimum pitch weight to be output}
    % MMAX = 12;		%{Maximum number of subharmonics in pitch calculation}
    % SPWEIGHT = 0.5;	%{Multiplier of spectral pitch weight}
    % maxfrequency = 5000; % maximum frequency for tone extraction
    
    WindowLength = round(fs/12.5);
    
    % index corresponding to max frequency for tone extraction
    maxfrequencyindex = ceil(1 + maxfrequency*WindowLength/fs);
    
    freq_all = fs*((1:round(WindowLength))'-1)/WindowLength;
    freq = freq_all(4:maxfrequencyindex);
    Offset = round(fs*timestep/1000);
    nwin = round((len-WindowLength)/Offset); % number of windows
    fftgain = 2^0.5 * rms(hann(WindowLength)) / WindowLength; % gain to be applied based on the FFT length
    t = ((0:nwin-1)' * Offset + 0.5 * WindowLength) ./ fs; % time vector for windows
    
    PitchPattern = zeros(nwin,MAXTONES);
    WeightPattern = zeros(nwin,MAXTONES);
    SpectralVirtual = true(nwin,MAXTONES);
    Salience = zeros(nwin,MAXTONES);
    PitchPatternScatter = zeros(nwin*MAXTONES,4);
    PitchPatternScatterIndex = 1;
    MIDISalience = zeros(nwin,128); % Quantized-pitch salience pattern using MIDI note numbers
    
    
    [PureTonalness, ComplexTonalness, Multiplicity, SumCompoundWeight,...
        SumSpectralWeight, SumVirtualWeight] = ...
            deal(zeros(nwin,1));
    
    % Process window loop
    for x = 1:nwin;
        
        
        
        start = round((x-1)*Offset + 1);
        finish = start + WindowLength - 1;
        Intensity_all = abs(fft(audio(start:finish)...
            * fftgain)).^2;
        L = pow2db(Intensity_all+1e-99);
        
        % The following peak identification and interpolation is the method
        % described in Terhardt et al.'s paper. We could introduce a better
        % approach here as an alternative.
        
        % Find peaks
        peakindex = (L(4:maxfrequencyindex) > L(3:maxfrequencyindex-1)) ...
            & (L(4:maxfrequencyindex) >= L(5:maxfrequencyindex+1)) ...
            & (L(4:maxfrequencyindex) - L(1:maxfrequencyindex-3) >= 7) ...
            & (L(4:maxfrequencyindex) - L(2:maxfrequencyindex-2) >= 7) ...
            & (L(4:maxfrequencyindex) - L(6:maxfrequencyindex+2) >= 7) ...
            & (L(4:maxfrequencyindex) - L(7:maxfrequencyindex+3) >= 7);
        
        % Interpolate to improve frequency estimate of peaks
        Intensity = Intensity_all(4:maxfrequencyindex);
        Lcrop = L(4:maxfrequencyindex);
        ToneF = freq(peakindex) + ...
            0.46*(Lcrop(circshift(peakindex,1)) ...
            -Lcrop(circshift(peakindex,-1)));
        ToneL = Lcrop(peakindex) + 1.6;
        NTones = length(ToneF);
        ToneRef = find(peakindex);
        
        
        % *****************************************
        % Calculate masking and pitch shift effects
        % (port from PsySound2 Pascal code)
        
        NTonesM = 0;
        toneBark = 13*atan(0.76*ToneF/1000) + 3.5*atan((ToneF/7500).^2);
        spectrumBark = 13*atan(0.76*freq/1000) + 3.5*atan((freq/7500).^2);
        % assume no more than 1000 tones...
        [LX, WS, FX, V] = deal(zeros(1000,1));
        for i = 1:NTones
            
            % Intensity of noise for each tone
            % paragraph after eq 7b
            Inoise = sum(Intensity(abs(toneBark(i)-spectrumBark)<=0.5 ...
                & (abs(ToneRef(i) - (1:length(spectrumBark))') > 2)));
            
            
            % the following should be vectorized!
            sumlo = 1e-99;
            sumhi = 1e-99;
            for j = 1:NTones
                if (j < i)
                    s = -24 - (0.23 / (ToneF(j)/1000)) + (0.2*ToneL(j)); %eq 7b
                    Lji = ToneL(j) - s * (spectrumBark(j) - toneBark(i));
                    sumlo = sumlo + 10.^(Lji/20);
                elseif (j > i)
                    Lji = ToneL(j) - 27 * (Fq2Bark(ToneF(j)) - toneBark(i));
                    sumhi = sumhi + 10.^(Lji/20);
                end; %if (ToneF [i] < ToneF [j]) then begin
            end; %for j := 1 to NTones do begin
            
            
            
            LXi = ToneL(i) - 10 * log10((sumlo + sumhi).^2 ...
                + Inoise + 10.^(Threshold(ToneF(i)) ./10)); %eq 4
            if LXi > 0
                NTonesM = NTonesM + 1;
                LX(NTonesM) = LXi;
                WS(NTonesM) = (1 - exp(LXi / -15)) ...
                    * (1 + 0.07 * (ToneF(i)/700 - 700/ToneF(i)).^2).^-0.5;
                FX(NTonesM) = ToneF(i);
                if PitchShift
                    LXid = ToneL(i) - 20*log10(sumlo);
                    LXidd = ToneL(i) - 20*log10(sumhi);
                    V(NTonesM) = 2e-4 * (ToneL(i) - 60) * ((ToneF(i)/1000) - 2) ...
                        + 1.5e-2 * exp(LXid/-20) * (3 - log(ToneF(i)/1000)) ...
                        + 3.0e-2 * exp(LXidd/-20) * (0.36 + log(ToneF(i)/1000)); %eq 10
                    FX(NTonesM) = ToneF(i) * (1 + V(NTonesM)); %eq 9
                end %if PitchShift
            end %if LXi > 0 then begin
        end %for i := 1 to NTones do begin
        
        if NTonesM > 0
            if NTonesM > 1
                % Sort pitches from high to low weight
                [WS,I] = sort(WS,'descend');
                LX = LX(I);
                FX = FX(I);
            end % if NTonesM > 1
            
            
            
            
            % Find virtual pitches
            delta = 0.08;
            MaxPitches = 0;
            while ~((MaxPitches+1 == NTonesM) || (WS(MaxPitches+1) < 0.7*WS(1)))
                MaxPitches = MaxPitches + 1;
            end
            
            [VPitch,VW] = deal(zeros(MaxPitches,MMAX));
            for i = 1:MaxPitches
                for m = 1:MMAX
                    Wim = 0;
                    for j = 1:NTonesM
                        if (j ~= i)
                            %n = trunc((m * FX(j))/ FX(i) + 0.5);
                            n = floor((m * FX(j))/ FX(i) + 0.5);
                            gamma = abs(n * FX(i) / (m * FX(j)) - 1);
                            if (gamma <= delta) && (n <= 20)
                                Cij = ((WS(i) * WS(j)) / (m * n)).^0.5 ...
                                    * (1 - (gamma / delta));
                                
                            else
                                Cij = 0;
                            end; %if (gamma <= delta) and (n <= 20) then begin
                            
                        else
                            Cij = 0;
                        end; %{if (j <> i then begin}
                        Wim = Wim + Cij;
                    end; %{for j := 1 to R do begin}
                    VPitch(i,m) = FX (i) / m;
                    if PitchShift
                        VPitch(i,m) = (VPitch(i,m)) * (1 + V(i) ...
                            - sign(m - 1) * 1e-3 * (18 + 2.5 * m ...
                            - (50 - 7 * m) * (VPitch(i,m) / 1000) ...
                            + 0.1 / (VPitch(i,m) / 1000)/(VPitch(i,m) / 1000)));
                        VW(i, m) = Wim / (1 + (((FX(i)/1000) / 0.8 / m).^4));
                    end
                end; %{for m:= 1 to MMAX do begin}
            end %{for i := 1 to MaxPitches do begin}
            
            isvirtual = 0; % whether or not there are any virtual pitches
            % Make compound pitch pattern
            CompoundPitch = zeros(MAXTONES,1);
            CompoundWeight = zeros(MAXTONES,1);
            Compoundsv = true(MAXTONES,1);
            PitchCount = 0;
            for i = 1:NTonesM
                if (WS(i) >= MINWEIGHT / SPWEIGHT)
                    PitchCount = PitchCount + 1;
                    CompoundPitch(PitchCount) = FX(i);
                    CompoundWeight(PitchCount) = WS(i) * SPWEIGHT;
                    PitchPatternScatter(PitchPatternScatterIndex,:) =...
                        [t(x), CompoundPitch(PitchCount), ...
                        CompoundWeight(PitchCount), 1];
                    PitchPatternScatterIndex = PitchPatternScatterIndex+1;
                end %{if WS [i] >= 2 *  MINWEIGHT then begin}
            end %{for i := 1 to NTonesM do begin}
            for i = 1:MaxPitches
                for m = 1:MMAX
                    if (VW(i, m) >= MINWEIGHT) && (VPitch(i, m) < 500)
                        isvirtual = 1;
                        PitchCount = PitchCount + 1;
                        CompoundPitch(PitchCount) = VPitch(i, m);
                        CompoundWeight(PitchCount) = VW(i, m);
                        Compoundsv(PitchCount) = false;
                        PitchPatternScatter(PitchPatternScatterIndex,:) =...
                            [t(x), CompoundPitch(PitchCount), ...
                            CompoundWeight(PitchCount), 0];
                        PitchPatternScatterIndex = PitchPatternScatterIndex+1;
                    end; %{if VW [i, m] >= MINWEIGHT then begin}
                end; %{for m := 1 to MMAX do begin}
            end; %{for i := 1 to MaxPitches do begin}
            
            if PitchCount > 0
                if PitchCount > 1
                    % Sort compound pitch pattern in terms of compound weight
                    [CompoundWeight,I] = sort(CompoundWeight,'descend');
                    CompoundPitch = CompoundPitch(I);
                    Compoundsv = Compoundsv(I);
                end %if PitchCount > 1
                
                if PitchCount > MAXTONES, PitchCount = MAXTONES; end
                PitchPattern(x,1:PitchCount) = CompoundPitch(1:PitchCount);
                WeightPattern(x,1:PitchCount) = CompoundWeight(1:PitchCount);
                SpectralVirtual(x,1:PitchCount) = Compoundsv(1:PitchCount);
                
                % Calculate Parncutt Measures
                Salience(x,1:PitchCount) = CompoundWeight(1:PitchCount) / ...
                    (CompoundWeight(1) * sum(CompoundWeight(1:PitchCount))).^0.5;
                PureTonalness(x) = (sum(CompoundWeight.^2)/5.2).^0.5;
                if isvirtual
                    ComplexTonalness(x) = max(CompoundWeight(~Compoundsv)/6.2);
                else
                    ComplexTonalness(x) = 0;
                end
                Multiplicity(x) = (sum(CompoundWeight)/CompoundWeight(1)).^0.5;
                
                SumCompoundWeight(x) = sum(CompoundWeight);
                SumSpectralWeight(x) = sum(CompoundWeight(Compoundsv));
                SumVirtualWeight(x) = sum(CompoundWeight(~Compoundsv));
                
                % Quantized pitch profile
                
                % MIDI = 12*log2(f/ATune) + 69;
                % Chroma = mod(MIDI-9,12);
                
                for i = 1:PitchCount
                    [MIDI,~,cents] = f2MIDI(CompoundPitch(i),1);
                    if MIDI > 0 && MIDI < 127
                        MIDISalience(x,MIDI+1) = MIDISalience(x,MIDI+1) + ...
                                Salience(x,i) * (100 - abs(cents))/100;
                            
                        if cents > 0
                            MIDISalience(x,MIDI+2) = MIDISalience(x,MIDI+2) + ...
                                Salience(x,i) * cents/100;
                            
                        elseif cents < 0
                            MIDISalience(x,MIDI) = MIDISalience(x,MIDI) + ...
                                Salience(x,i) * -cents/100;
                        end
                    end
                end
            end % PitchCount > 0
        end
        
        %disp('hello')
        % *****************************************
        
    end % window loop
    
    % Salience of chroma
    
    % Time-averaged MIDI salience
    MIDISalienceMean = mean(MIDISalience);
    
    ChromaSalience = zeros(nwin,12);
    for c = 1:12
        ChromaSalience(:,c) = ...
            sum(MIDISalience(:,c+[9,21,33,45,57,69,81,93,105]),2);
        % 'A' 'A#' 'B' 'C' 'C#' 'D' 'D#' 'E' 'F' 'F#' 'G' 'G#'
    end
   
    if chromaorder
         % change to cycle-of-fifths order
        ChromaSalienceTemp = ChromaSalience;
        ChromaSalience(:,2) = ChromaSalienceTemp(:,8);
        ChromaSalience(:,4) = ChromaSalienceTemp(:,10);
        ChromaSalience(:,6) = ChromaSalienceTemp(:,12);
        ChromaSalience(:,8) = ChromaSalienceTemp(:,2);
        ChromaSalience(:,10) = ChromaSalienceTemp(:,4);
        ChromaSalience(:,12) = ChromaSalienceTemp(:,6);
        chomalabel = {'A' 'E' 'B' 'F#' 'C#' 'G#' 'D#' 'A#' 'F' 'C' 'G' 'D'};
    else
        chromalabel = {'A' 'A#' 'B' 'C' 'C#' 'D' 'D#' 'E' 'F' 'F#' 'G' 'G#'};
    end
    
    % Time-averaged chroma salience
    ChromaSalienceMean = mean(ChromaSalience);
    ChromaStats = zeros(18,12);
    for n = 1:12
        ChromaStats(:,n) = pitchstats(ChromaSalience(:,n));
    end
    
    
    % MIDI salience plots
    figure('Name', 'MIDI note number salience')
    subplot(2,1,1)
    imagesc(t,1:128,MIDISalience')
    xlabel('Time (s)');
    ylabel('MIDI note number')
    set(gca,'TickDir', 'out');
    set(gca,'YDir','normal')
    ylim([1 128]);
    
    subplot(2,1,2)
    bar(1:128,MIDISalienceMean)
    xlabel('MIDI note number')
    xlim([1 128])
    
    
    % Chroma salience plots
    figure('Name','Chroma salience')
    subplot(2,1,1)
    imagesc(t,1:12,ChromaSalience')
    xlabel('Time (s)');
    ylabel('Chroma')
    set(gca,'YTick',1:12)
    set(gca,'YTickLabel',chomalabel)
    set(gca,'TickDir', 'out');
    set(gca,'YDir','normal')
    %legend('show');
    ylim([1 12]);
    
    subplot(2,1,2)
    bar(1:12,ChromaSalienceMean)
    xlabel('Chroma number')
    %xlim([1 12])
    set(gca,'XTickLabel',chomalabel)
    
    
    
    % Compare pitch profile with chord and key templates
    order = 6; % this is the amount of smoothing applied in chord change estimation
    %[chord2,cLH,ChrdChTimes,chroTSpec] = PitchProfile(Salience,timePoints,order)
    
    
    
    % Pitch 3-D scattergram
    PitchPatternScatter = PitchPatternScatter(1:PitchPatternScatterIndex,:);
    
    VP = find(PitchPatternScatter(:,4)==0);
    SP = find(PitchPatternScatter(:,4)==1);
    figure('Name','Pitch Pattern');
    if ~isempty(SP)
        scatter3(PitchPatternScatter(SP,1), ...
            PitchPatternScatter(SP,2), ...
            PitchPatternScatter(SP,3), ...
            'b')
        hold on
    end
    if ~isempty(VP)
        scatter3(PitchPatternScatter(VP,1), ...
            PitchPatternScatter(VP,2), ...
            PitchPatternScatter(VP,3), ...
            'r')
    end
    
    % Pitch time-series metrics
    figure('Name', 'Pitch Time-Series Metrics')
    subplot(3,1,1)
    plot(t,PureTonalness,'r')
    xlabel('Time (s)');
    ylabel('Pure Tonalness')
    
    subplot(3,1,2)
    plot(t,ComplexTonalness,'r')
    xlabel('Time (s)');
    ylabel('Complex Tonalness')
    
    subplot(3,1,3)
    plot(t,Multiplicity,'r')
    xlabel('Time (s)');
    ylabel('Multiplicity')
    
    
    % Tables of means
    fig1 = figure('Name','Virtual Pitch Statistics over Time');
    
    % Pitch parameter stats
    [PureTonalnessStats,labels] = pitchstats(PureTonalness);
    ComplexTonalnessStats = pitchstats(ComplexTonalness);
    MultiplicityStats = pitchstats(Multiplicity);
    PeakSalienceStats = pitchstats(max(Salience,[],2));
    SumSalienceStats = pitchstats(sum(Salience,2));
    SumPitchWeightStats = pitchstats(SumCompoundWeight);
    SumSpectralWeightStats = pitchstats(SumSpectralWeight);
    SumVirtualWeightStats = pitchstats(SumVirtualWeight);
    
    data = [PureTonalnessStats ComplexTonalnessStats MultiplicityStats PeakSalienceStats SumSalienceStats SumPitchWeightStats SumSpectralWeightStats SumVirtualWeightStats];
    
    table1 =  uitable('Data',data,...
                      'RowName',labels,...
                      'ColumnName',{'Pure Tonalness' 'Complex Tonalness' 'Multiplicity' 'Highest Salience' 'Salience Sum' 'Compound Weight Sum' 'Spectral Weight Sum' 'Virtual Weight Sum'});
    
    
    table2 = uitable('Data',ChromaStats,...
                     'RowName',labels,...
                     'ColumnName',chomalabel);
    
    [~,table] = disptables(fig1,[table1 table2]);
    OUT.tables = table;
    
    %disp('hello')
    % *** CREATING A TABLE OR MULTIPLE TABLES ***
    % You may include tables to display your results using AARAE's
    % disptables.m function - this is just an easy way to display the
    % built-in uitable function in MATLAB. It has several advantages,
    % including:
    %   * automatically sizing the table(s);
    %   * allowing multiple tables to be concatenated;
    %   * allowing concatenated tables to be copied to the clipboard
    %     by clicking on the grey space between the tables in the figure;
    %   * if its output is used as described below, returning data to
    %     the AARAE GUI in a format that can be browsed using a bar-plots.
    %   * and values from tables created this way are written into the log
    %     file for the AARAE session, which provides another way of
    %     accessing the results together with a complete record of the
    %     processing that led to the results.
    
    %         fig1 = figure('Name','Virtual Pitch Analysis');
    %         table1 = uitable('Data',[t PureTonalness ComplexTonalness Multiplicity Salience],...
    %                     'ColumnName',{'Time (s)','Pure Tonalness','Complex Tonalness','Multiplicity','Salience'},...
    %                     'RowName',{'Results'});
    %     disptables(fig1,table1); % see below for a better alternative!
    % It is usually a very good idea to return your table(s) data to AARAE
    % as a function output field. Doing this creates a table leaf in the
    % Results section of the GUI, which can be explored as bar plots within
    % the GUI, and also writes the table contents to the log file (as a
    % comma delimited table that can easily be interpreted by a
    % spreadsheet). To do this, simply use the output of the
    % disptables funtion as follows:
    %     [~,table] = disptables(fig1,table1);
    % And include your table in the output data structure
    %      OUT.tables = table;
    %
    %
    % If you have multiple tables to combine in the figure, you can
    % concatenate them:
    %       disptables(fig1,[table3 table2 table1]);
    % (Note that the concatenation is in descending order - i.e. the first
    % listed table in the square brackets will be displayed at the bottom,
    % and the last listed table will be displayed at the top.)
    %
    % You may export these tables to be displayed as bar plots as if you
    % were doing it for a single table:
    %       [~,tables] = disptables(fig1,[table3 table2 table1]);
    %       OUT.tables = tables;
    %
    % The disptables function will take care of allocating each table to a
    % different barplot, there is no need to generate more than one .tables
    % field to display all your tables.
    
    
    
    % *** MAKING PLOTS ***
    % You may also include figures to display your results as plots.
    %     t = linspace(0,duration,length(audio));
    %     figure;
    %     plot(t,audio);
    % All figures created by your function are stored in the AARAE
    % environment under the results box.
    
    % disp('hello')
    
    % *** CREATING A RESULTS LEAF (FOR BIG NON-AUDIO DATA) ***
    % You may output data to be plotted in a variety of charts, including
    % lines, mesh, surf, imagesc, loglog, and others depending on the
    % number of dimensions of your data using the doresultleaf.m function:
    % E.g.:
    %
    %         doresultleaf(myresultvariable,'Type [units]',{'Time','channels','Frequency'},...
    %                      'Time',      t,                's',           true,...
    %                      'channels',  chanID,           'categorical', [],...
    %                      'Frequency', num2cell(bandfc), 'Hz',          false,...
    %                      'name','my_results');
    %     %
    % Input arguments:
    % #1: Your data variable. It can be multidimensional, make sure you
    %     specify what each dimension is.
    % #2: What is your data variable representing? is it level? is it
    %     reverb time? make sure you label it appropriately and assign
    %     units to it, this second argument is a single string.
    % #3: This is a cell array where each cell contains the name of each
    %     dimension.
    %
    % #4: Name of your first dimension. (String)
    % #5: Matrix that matches your first dimension, in this case time.
    % #6: Units of your first dimension. (String)
    % #7: Can this dimension be used as a category? (true, false, [])
    %
    % Replicate arguments 4 to 7 for as many dimensions as your data
    % variable has.
    %
    % The second last input argument is the string 'name', this helps the
    % function identify that the last input argument is the name that will
    % be displayed in AARAEs categorical tree under the results leaf.
    
    
    
    % And once you have your result, you should set it up in an output form
    % that AARAE can understand.
    if isstruct(IN)
        
        % *** OUTPUTTING AUDIO ***
        % Most analysers do not output audio. If you wish to output audio,
        % first consider whether your analyser should be a processor
        % instead. If it should be an analyser, then audio can be output as
        % follows:
        %OUT = IN; % You can replicate the input structure for your output
        %OUT.audio = audio; % And modify the fields you processed
        % However, for an analyser, you might not wish to output audio (in
        % which case the two lines above might not be wanted.
        
        
        
        % *** OUTPUTTING NEW FIELDS ***
        % Or simply output the fields you consider necessary after
        % processing the input audio data, AARAE will figure out what has
        % changed and complete the structure. But remember, it HAS TO BE a
        % structure if you're returning more than one field:
        
        %         OUT.properties.duration = duration;
        %         OUT.properties.maximum = maximum;
        %         OUT.properties.minimum = minimum;
        % The advantages of providing outputs as subfields of properties
        % are that AARAE has a button that opens a window to display
        % properties, and that properties are also written to the log file
        % of the AARAE session. Outputing fields without making them
        % subfields of properties is possible, but this makes the output
        % data harder to access.
        % Note that the above outputs might be considered to be redundant
        % if OUT.tables is used, as described above. Generally the
        % properties fields output is suitable for small data only (e.g.
        % single values). Tables is best for small and medium data, while a
        % results leaf is is best for big data.
        
        
        
        % *** FUNCTION CALLBACKS ***
        % The following outputs are required so that AARAE can repeat the
        % analysis without user interaction (e.g. for batch processing),
        % and for writing the log file.
        OUT.funcallback.name = 'Terhardt_VirtualPitch.m'; % Provide AARAE
        % with the name of your function
        OUT.funcallback.inarg = {timestep,PitchShift,MINWEIGHT,MMAX,SPWEIGHT,maxfrequency,ATune,chromaorder,cal,fs};
        % assign all of the input parameters that could be used to call the
        % function without dialog box to the output field param (as a cell
        % array) in order to allow batch analysing. Do not include the
        % first (audio) input here.
        
        
        
        % AARAE will only display the output in the main GUI if it includes
        % tables, audio or other fields. It will not display if only
        % function callbacks are returned. Result leaves are not created by
        % the function's main output, but instead are created by calling
        % AARAE's doresultleaf as described above, and these will be
        % displayed in AARAE's main GUI regardless of the function's main
        % outputs.
    else
        % You may increase the functionality of your code by allowing the
        % output to be used as standalone and returning individual
        % arguments instead of a structure.
        %OUT = audio;
    end
    
else
    % AARAE requires that in case that the user doesn't input enough
    % arguments to generate audio to output an empty variable.
    OUT = [];
end
end % eof

function L = Threshold(f)
% hearing threshold
% input f is frequency in Hz
% output L is threshold in dB
f=f/1000;
L = 3.64 * f.^-0.8 ...
    - 6.5 * exp(-0.6 * (f - 3.3).^2) ...
    + 1e-3 * f.^4;
end

function B = Fq2Bark(f)
% critical band rate corresponding to a given frequency
% input f is frequency in Hz
% output B is critical band rate in Barks
f=f/1000;
B = 13 * atan(0.76 * f) + 3.5 * atan ((f/7.5).^2);
end


function[MIDI, chroma, cents] = f2MIDI(f,method,ATune)
% calculates the MIDI note number and chroma corresponding to an input
% frequency
% ATune is the tuning frequency of A (usually 440 Hz)

if ~exist('ATune','var'), ATune = 440; end
if ~exist('method','var'), method = 1; end

switch method
    
    case 0 % output unrounded values
        MIDI = 12*log2(f/ATune) + 69;
        chroma = mod(MIDI-10,12)+1;
        cents = 100*(MIDI-round(MIDI));
        
    case 1 % output rounded vaues
        MIDI = 12*log2(f/ATune) + 69;
        cents = 100*(MIDI-round(MIDI));
        MIDI = round(MIDI);
        chroma = mod(MIDI-10,12)+1; 
        
end
end
                

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pitch Quantization Functions
%
function out = createSpectrum(pitch, notes)
% get the frequencies in this frame - convert to midi
notesInFrame = log2(pitch(:,1));

% get the corresponding saliences
saliencesInFrame = pitch(:,2);

% setup vector
out = zeros(1, length(notes));
for j =1:length(notesInFrame)
    % which two are we looking at
    lowNoteIndex  = find(notes < notesInFrame(j), 1, 'last');
    highNoteIndex = find(notes > notesInFrame(j), 1, 'first');
    
    % what's the difference
    difference = notes(highNoteIndex) - notes(lowNoteIndex);
    
    % Multiply by distance from each of the notes
    out(lowNoteIndex)  = saliencesInFrame(j) * ...
        ((notes(highNoteIndex) - notesInFrame(j))/difference);
    
    out(highNoteIndex) = saliencesInFrame(j) * ...
        ((notesInFrame(j) - notes(lowNoteIndex))/difference);
end
end %createSpectrum





function out = createChromaPattern(notespectrum)
% output with two columns: chroma, sum of salience
chromaoffset = 0; % integer to change numerical assignment of chroma
out = zeros(1,12); % set up vector
chromapattern = zeros(2,12);
% chromapattern(:,1) = 1:12 - 1 + chromaoffset; % chroma in simple assending order
chromapattern(1:12,1) = mod((1:12 - 1 + chromaoffset) * 7, 12); % chroma in cycle of fifths order
% sum the saliences of each chroma
for j = 1:length(notespectrum)
    chromapattern((mod(j,12)+1),2) = chromapattern((mod(j,12)+1),2)...
        + notespectrum(j);
    chromapattern = sortrows(chromapattern,1);
    out(1:12,1) = chromapattern(1:12,2);
end %end createChromaPattern
end



function [data,meaning] = pitchstats(in)
meaning = {'mean','standard deviation','maximum','99%','98%','97%','96%','95%',...
    '90%','80%','70%','60%','median','40%','30%','20%','10%','minimum'};
data = [mean(in);...
        std(in);...
    max(in);...
    prctile(in,99);...
    prctile(in,98);...
    prctile(in,97);...
    prctile(in,96);...
    prctile(in,95);...
    prctile(in,90);...
    prctile(in,80);...
    prctile(in,70);...
    prctile(in,60);...
    median(in);...
    prctile(in,40);...
    prctile(in,30);...
    prctile(in,20);...
    prctile(in,10);...
    min(in)];

end


%**************************************************************************
% Copyright (c) 2015, Densil Cabrera
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