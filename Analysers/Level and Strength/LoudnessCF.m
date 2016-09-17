function OUT = LoudnessCF(IN,HL,ov,fs,cal)
% This function calculates loudness time-varying loudness, specific
% loudness, loudness fluctuation, and sharpness using the loudness model of
% Chalupper and Fastl:
% J. Chalupper & H. Fastl (2002) "Dynamic Loudness Model (DLM) for Normal
% and Hearing-Impaired Listeners" Acta Acustica United with Acustica,
% 88:378-386, 2002.
% J. Chalupper (2000) "Modellierung der Lautstärkeschwankung für Normalund
% Schwerhörige" in DAGA 2000, 254-255.
%
% INPUT ARGUMENTS
% IN - audio signal (1 channel only)
% HL - hearing loss in dB for 24 critical bands (Bark)
% ov - time window overlap in percentage
% When run directly from the AARAE GUI, HL and ov are zeros. If you wish to
% use non-zero inputs while using the AARAE GUI, one way is to write a
% workflow function (by adapting code generated in AARAE's log file)
% fs and cal inputs are only used if a non-structure (vector or matrix) is
% the primary input.
%
% OUTPUTS
% * Time-varying Loudness
% * Time-varying Sharpness (Zwicker & Fastl)
% * Time-averaged Specific Loudness
% * Time-averaged Loudness statistics, including fluctuation rate
% * Time-averaged Sharpness statistics

% Author: Josef Chalupper (josef.chalupper@siemens.com) original version: 12.12.2000
% New version (with comments and examples by the PsySound team): 6.1.2007 for psysound3
% Adjusted version by Ella Manor for AARAE 04-09-2015

% *************************************************************************

% *************************************************************************
if isstruct(IN)
    
    % only 1 channel analysis at present:
    IN = choose_from_higher_dimensions(IN,1,1);
    % or the following for multichannel, once it is implemented
    %IN = choose_from_higher_dimensions(IN,2,1); 

    
    
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    if isfield(IN,'cal') % Get the calibration offset if it exists
        cal = IN.cal;
    else
        warndlg('Calibration data missing - please calibrate now.','LoudnessCF','modal');
        IN = cal_aarae(IN);
        if isempty(IN) % user pressed cancel
            OUT = IN;
            return
        end
        cal = IN.cal;
        %         h=warndlg('Calibration data missing - please calibrate prior to calling this function.','LoudnessCF','modal');
        %         uiwait(h)
        %         OUT = []; % you need to return an empty output
        %         return % get out of here!
    end
    
    % not used until multichannel is implemented
    if isfield(IN,'chanID') % Get the channel ID if it exists
        chanID = IN.chanID;
    else
        % or make a chanID using AARAE's utility function
        chanID = makechanID(size(audio,2),0);
    end
    
    if isfield(IN,'name') % Get the AARAE name if it exists
        name = IN.name;
    else
        name = '';
    end
else
    audio = IN;
    name = '';
    chanID = makechanID(size(audio,2),0);
end
% *************************************************************************

if ~isempty(audio) && ~isempty(fs) && ~isempty(cal)
    [~, nchan, bands] = size(audio);
    if bands > 1
        audio = sum(audio,3);
        disp('Multiband audio has been mixed.')
    end
    
    if nchan > 1
        audio = audio(:,1:2);
        if nchan > 2
            disp('Only the first two channels are analysed by LoudnessCF.')
            nchan = 2;
        end
    end
    
    % signal calibration offset
    calconstant = 109.8999; % by validating 1 kHz at 40 dB = 1 sone
    cal = cal-calconstant;
    if length(cal) == 1
        audio = audio .* 10.^(cal/20);
    else
        audio(:,1) = audio(:,1) .* 10.^(cal(1)/20);
        audio(:,2) = audio(:,2) .* 10.^(cal(2)/20);
    end
    
    %disp(['rms level of the entire wave ', num2str(10*log10(mean(audio.^2)+10e-99)+calconstant), ' dB'])
    
    % *************************************************************************
    % MAIN analyser starts here
    % *************************************************************************
    
    
    % ********* set up the analysis parameters etc ************************
    f_abt = 1/2e-3; % 2 ms sampling period
    
    % window size
    [t_pa,w,t_sb,t_sa,t_pb] = staticParamDLM;
    [h,~, ~] = tep_window(t_pb,t_pa,t_sb,t_sa,w,fs);
    wl = length(h);
    if mod(wl, 2)
        wl = wl + 1;
    end
    
    % overlap
    if ~exist('ov','var')
        ov = 0;
    elseif isempty(ov),
        ov = 0;
    end
    
    ov_samp = wl*ov/100;
    offset = wl - ov_samp;
    period = offset ./ fs;
    wr = round(1 / period,2);
    
    % Hearing loss
    if ~exist('HL','var')
        HL = zeros(1,24);
    elseif isempty(HL),
        HL = zeros(1,24);
    end
    
    k  = 0.8; % HL factor
    % Splitting of hearing loss into HL_ihc and HL_ohc
    HL_ohc = k.*HL;
    HL_ihc = HL-HL_ohc;
    
    % Approximation of the transfer function through the human outer
    % and and middle ear
    [b, a] = butter_hp(fs); % generate the butterworth filter coeffs
    
    % filter state vector.
    Z = [];
    
    % Calculation of coefficients of critical band filterbank
    S = make_fttbank1(fs);
    
    kern_l = [];
    
    % Smoothed critical band loudness filter creation
    [smooth.b, smooth.a] = int_tp(f_abt);
    smooth.Zfa = [];
    smooth.Zfb = [];
    
    % *************** analyse the audio ***********************************
    
    
    % Run the butterworth filter
    [sig, ~] = filter(b, a, audio, Z);
    
    % Applying critical band filterbank
    [fgrp, ~] = ftt_bank1(sig, S, f_abt,fs);
    fgrp_d    = damp_a0(fgrp, HL_ihc); % Attenuation due to outer & middle
    % ear and inner hair cell hearing
    % loss
    
    % Calculation of main loudness
    kern_l = [kern_l; kernlaut24_two(fgrp_d, HL_ohc)];
    
    % Calculation of forward masking (aka "post masking")
    try
        kern_dyn = post_maskn(kern_l, f_abt); % no effect = kern_l
    catch
        % caught the case where no postprocessing is possible, use a string to
        % indicate to the calling function.
        N = 'no postprocessing';
        main_N = 0;
        spec_N = 0;
    end
    
    % Calculation of spectral masking and spectral summation of
    % specific loudness labelled N'
    [spec_N, lauth] = flankenlautheit24(kern_dyn);
    
    % summation of critical band loudness
    kl = bark_sum(spec_N); % The integral of specific loudness over critical-band rate
    
    % Smoothed critical band loudness (specific loudness)
    [main_N, smooth.Zfa] = filter(smooth.b, smooth.a, kl, smooth.Zfa);
    %main_N(find(main_N < 0)) = 0;
    main_N(main_N < 0) = 0;
    
    % Loudness integration
    [N, smooth.Zfb] = filter(smooth.b, smooth.a, lauth, smooth.Zfb);
    %N(find(N < 0)) = 0;
    N(N<0)=0;
    
    % calculate sharpness - over time
    [r, ~] = size(spec_N);
    acum = zeros(r, 1);
    for i = 1:r
        acum(i) = sharpness_Fastl(spec_N(i,:));
    end
    
    % calculate overall fluctuation
    % Note that this is loudness fluctuation (as described by Chalupper
    % & Fastl, NOT fluctuation strenth (as per Zwicker & Fastl) - to
    % answer one of the most commonly asked questions since the
    % PsySound3 implementation of this.
    lf=fluct(main_N);
    
    % *************************************************************************
    % Prep data for visualisation
    % *************************************************************************
    
    % Nominal
    tPeriod = 2e-3; % 2 ms
    timePoints = (0:length(N)-1)' * tPeriod;
    
    % Time-averaged specific loudness
    T_spec_N = mean(spec_N);
    
    OUT=[];
    
    % *************************************************************************
    % Data Presentation
    % *************************************************************************
    
    % ********* TABLES *********
    
    % Loudness statistics, adopted from Loudness_MG2b code
    Nmean = mean(N);
    Nstd = std(N);
    Nmax = max(N);
    N1 = prctile(N,99);
    N2 = prctile(N,98);
    N3 = prctile(N,97);
    N4 = prctile(N,96);
    N5 = prctile(N,95);
    N10 = prctile(N,90);
    N20 = prctile(N,80);
    N30 = prctile(N,70);
    N40 = prctile(N,60);
    N50 = median(N);
    N60 = prctile(N,40);
    N70 = prctile(N,30);
    N80 = prctile(N,20);
    N90 = prctile(N,10);
    Nmin = min(N);
    
    dataN = [Nmean;Nstd;Nmax;lf;N1;N2;N3;N4;N5;N10;N20;N30;N40;N50;N60;N70;N80;N90;Nmin];
    
    % Sharpness statistics, adopted from Loudness_MG2b code
    Smean = mean(acum);
    Sstd = std(acum);
    Smax = max(acum);
    S1 = prctile(acum,99);
    S2 = prctile(acum,98);
    S3 = prctile(acum,97);
    S4 = prctile(acum,96);
    S5 = prctile(acum,95);
    S10 = prctile(acum,90);
    S20 = prctile(acum,80);
    S30 = prctile(acum,70);
    S40 = prctile(acum,60);
    S50 = median(acum);
    S60 = prctile(acum,40);
    S70 = prctile(acum,30);
    S80 = prctile(acum,20);
    S90 = prctile(acum,10);
    Smin = min(acum);
    
    dataS = [Smean;Sstd;Smax;S1;S2;S3;S4;S5;S10;S20;S30;S40;S50;S60;S70;S80;S90;Smin];
    
    % generate tables of results
    
    fig1 = figure('Name','Time-varying Dynamic Loudness (C&F) and Sharpness Statistics');
    table1 = uitable('Data',dataN,...
        'ColumnName',{'Loudness'},...
        'RowName',{'Mean','Standard deviation','Maximum','Fluctuation',...
        'N1','N2','N3','N4',...
        'N5','N10','N20','N30','N40','N50 (median)','N60',...
        'N70','N80','N90','Minimum'});
    table2 = uitable('Data',dataS,...
        'ColumnName',{'Sharpness'},...
        'RowName',{'Mean','Standard deviation','Maximum',...
        'S1','S2','S3','S4',...
        'S5','S10','S20','S30','S40','S50 (median)','S60',...
        'S70','S80','S90','Minimum'});
    
    [~,tables] = disptables(fig1,[table1 table2]); % AARAE function
    
    
    
    OUT.tables = tables;
    
    % ********* CHARTS *********
    
    y_lim = max(acum)+1;
    spec_vec = 1:240;
    % Figure for charts
    figure('Name',['Loudness (C&F) of ',name])
    
    subplot(2,2,1:2)
    
    % Time-varying loudness and sharpness
    [ax,line1,line2] = plotyy(timePoints,N,timePoints,acum);
    title ('Time-Varying Loudness and Sharpness');
    xlabel('Time (s)');
    ax(1).YColor = 'k';
    ax(2).YColor = 'r';
    ylabel(ax(1), 'Loudness (sone)','Color','k');
    ylabel(ax(2), 'Sharpness (acum)','Color','r');
    set(ax(2),'YLim',[0 y_lim]);
    set(line1,'Color','k','DisplayName','loudness');
    set(line2,'Color','r','DisplayName','sharpness');
    legend('show','Location','NorthEast');
    
    % Time-averaged specific loudness as a fucntion of critical band
    % figure
    subplot(2,2,4)
    plot(1:240,T_spec_N,'k-');
    ax=gca;
    ax.Title.String = 'Time-Averaged Specific Loudness';
    ax.XLabel.String = 'Critical Band Rate (Bark)';
    ax.XLim = [0 length(T_spec_N)+10];
    ax.XTickLabel = {'0','5','10','15','20','25'};
    ax.YLabel.String = 'Loudness (sones/Bark)';
    hold off;
    
    % Specific loudness spectrogram
    subplot(2,2,3)
    imagesc(timePoints, 0.1:0.1:24, spec_N');
    % cH = colorbar;
    set(gca,'YDir','normal');
    ax=gca;
    axis tight;
    ax.Title.String = 'Specific Loudness';
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Critical Band Rate (Bark)';
    hold off;
    
    
    % ******** AARAE RESULTS LEAVES **********
    if isstruct(IN)
        % Time-varying loudness results leaf
        doresultleaf(N,'Loudness [sone]',{'time'},...
            'Time',timePoints','s',true,...
            'loudnesstype', {'Loudness over time'}, 'categorical',[],...
            'name','Time_varying_loudness');
        
        % Time-varying sharpness results leaf
        doresultleaf(acum,'Sharpness [acum]',{'time'},...
            'Time',timePoints','s',true,...
            'loudnesstype', {'Sharpness over time'}, 'categorical', [],...
            'name','Time_varying_sharpness');
        
        % Time-averaged Specific loudness results leaf
        doresultleaf(T_spec_N','Specific Loudness [sones/Bark]',{'Critical Band Rate'},...
            'Critical Band Rate',[0.1:0.1:24]','Bark',true,...
            'loudnesstype', {'Loudness over critical band'}, 'categorical', [],...
            'name','Time_averaged_specific_loudness');
        
        % Specific loudness spectrogram results leaf
        doresultleaf(spec_N','Specific Loudness [sones/Bark]',{'time','Critical Band Rate'},...
            'Critical Band Rate',[0.1:0.1:24]','Bark',true,...
            'Time',timePoints','s',true,...
            'loudnesstype', {'Loudness over critical band'}, 'categorical', [],...
            'name','Time_varying_specific_loudness');
        
    end
    
    OUT.funcallback.name = 'LoudnessCF.m';
    OUT.funcallback.inarg = {HL,ov,fs,cal};
    
else
    OUT = [];
end


%**************************************************************************
% wrapper Copyright (c) 2015, Ella Manor & Densil Cabrera
% Main code used with permission from Joseph Chalupper
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