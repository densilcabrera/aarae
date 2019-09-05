function out = LateralFractions(data,fs,startthresh,threshmethod,bpo,omniweight,doplot)
% This function calculates spatial fractions based on first order Ambisonic
% mic IRs measured in an auditorium.
%
% Lateral fractions are calculated as per ISO 3382-1, but fractions are
% also calculated for the two other first-order channels (front-back and
% up-down), using the same time periods.
%
% Inputs:
% data is the audio structure or waveform
%
% fs is sampling rate in Hz
%
% startthresh is the threshold in dB used to detect the truncation point
% just before the direct sound impulse
%
% bpo is bands per octave (1 for octave band, 3 for 1/3-octave band)
%
% omniweight is the relative weight of the zeroth vs 1st order channels
% (e.g. chan 1 may have -2.4 dB sensitivity relative to on-axis chans 2-4).
% This depends on the Ambisonics encoding format.
%
% doplot is used to choose whether a chart is generated or not



if nargin < 6, doplot = 1; end
if nargin < 5, omniweight = -2.4; end
if nargin < 4, bpo = 1; end
if nargin < 3
    startthresh = -20;
    %dialog box for settings
    prompt = {'Threshold for IR start detection', ...
        'Detect start before filtering (0) or after filtering (1)',...
        'Bands per octave (1 | 3)', ...
        'Weight of omni vs fig-8 channels (dB)',...
        'Plot (0|1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {num2str(startthresh),'0','1','-2.4','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        startthresh = str2num(answer{1,1});
        threshmethod = str2num(answer{2,1});
        bpo = str2num(answer{3,1});
        omniweight = str2num(answer{4,1});
        doplot = str2num(answer{5,1});
    else
        out = [];
        return
    end
end
if isstruct(data)
    data = choose_from_higher_dimensions(data,3,1);
    IR = data.audio;
    fs = data.fs;
else
    IR = data;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
            'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

if ~isempty(IR) && ~isempty(fs) && ~isempty(startthresh) && ~isempty(bpo) && ~isempty(doplot)
    
    %--------------------------------------------------------------------------
    % TRUNCATION PRIOR TO FILTERING
    %--------------------------------------------------------------------------
    
    % Check the input data dimensions
    [len, chans, bands] = size(IR);
    if bands>1
        warndlg('The function is not suitable for multiband audio. Multiband audio has been mixed down, but results should be interpreted with caution.')
        IR = mean(IR,3);
    end
    
    if chans > 4
        IR = IR(:,1:4);
        %warndlg('Only the first four channels are analysed, assuming 1st order Ambisonics format')
        chans = 4;
    end
    if chans < 4
        warndlg('Currently this analyser requires 4 channel input in 1st order Ambisonics format')
        out = [];
        return
    end
    
    if len/fs < 1
        warndlg('The impulse response is less than 1 s. Zero padding has been applied, but results should be interpreted with caution.')
        IR = [IR; zeros(round(fs),chans)];
    end
    
    % Adjust IR for equal channel sensitivity
    IR(:,1) = IR(:,1) .* db2mag(-omniweight);
    
    if threshmethod == 0
        % auto-crop start of individual audio columns
        IR = autocropstart_aarae(IR,startthresh,1);
        
        % Truncate sections of the IR prior to filtering
        Direct10 = IR(1:1+floor(fs*0.01),:); % Truncate first 10 ms
        Early5_80 = IR(1+floor(fs*0.005):1+floor(fs*0.08),:);
        Early0_80 = IR(1:1+floor(fs*0.08),:);
        Late80_end = IR(1+floor(fs*0.08):end,:);
        
        %         Direct10 = IR(1:1+floor(fs*0.01),:); % Truncate first 10 ms
        %         Early20_100 = IR(1+floor(fs*0.02):1+floor(fs*0.1),:); % 20-100 ms
        %         Late100_200 = IR(1+floor(fs*0.1):1+floor(fs*0.2),:); % 100-200 ms
        %         All20_1000 = IR(1+floor(fs*0.02):1+floor(fs*1),:); % 20-1000 ms
        %         Late100_1000 = IR(1+floor(fs*0.1):1+floor(fs*1),:); % 100-1000 ms
        
        
        %--------------------------------------------------------------------------
        % FILTERING
        %--------------------------------------------------------------------------
        
        % Use inbuilt AARAE filterbank
        if bpo == 3
            bandfc = [100,125,160,200,250,315,400,500,630,800,1000,1250,1600, ...
                2000,2500,3150,4000,5000];
            All_IR = thirdoctbandfilter_viaFFT(IR,fs,bandfc,15,round(0.1*fs)); % IR
            Direct10 = thirdoctbandfilter_viaFFT(Direct10,fs,bandfc,15,round(0.1*fs));
            Early5_80 = thirdoctbandfilter_viaFFT(Early5_80,fs,bandfc,15,round(0.1*fs));
            Early0_80 = thirdoctbandfilter_viaFFT(Early0_80,fs,bandfc,15,round(0.1*fs));
            Late80_end = thirdoctbandfilter_viaFFT(Late80_end,fs,bandfc,15,round(0.1*fs));
            
            
        else
            bandfc = [125,250,500,1000,2000,4000];
            All_IR = octbandfilter_viaFFT(IR,fs,bandfc,12,round(0.1*fs)); % IR
            Direct10 = octbandfilter_viaFFT(Direct10,fs,bandfc,12,round(0.1*fs));
            Early5_80 = octbandfilter_viaFFT(Early5_80,fs,bandfc,12,round(0.1*fs));
            Early0_80 = octbandfilter_viaFFT(Early0_80,fs,bandfc,12,round(0.1*fs));
            Late80_end = octbandfilter_viaFFT(Late80_end,fs,bandfc,12,round(0.1*fs));
            
            
        end
    end
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    % TRUNCATION AFTER FILTERING
    %--------------------------------------------------------------------------
    
    if threshmethod == 1
        % apply filterbank
        if bpo == 3
            bandfc = [100,125,160,200,250,315,400,500,630,800,1000,1250,1600, ...
                2000,2500,3150,4000,5000];
            IR = thirdoctbandfilter_viaFFT(IR,fs,bandfc,15,round(0.1*fs)); % IR
            
        else
            bandfc = [125,250,500,1000,2000,4000];
            IR = octbandfilter_viaFFT(IR,fs,bandfc,12,round(0.1*fs)); % IR
        end
        % autocrop start of individual audio columns (removes band time
        % alignment)
        IR = autocropstart_aarae(IR,startthresh,1);
        
        % Truncate sections of the IR
        Direct10 = IR(1:1+floor(fs*0.01),:,:); % Truncate first 10 ms
        Early5_80 = IR(1+floor(fs*0.005):1+floor(fs*0.08),:,:);
        Early0_80 = IR(1:1+floor(fs*0.08),:,:); % 0-80 ms
        Late80_end = IR(1+floor(fs*0.08):end,:,:); % 80-end ms
        All_IR = IR;
    end
    
    
    
    %--------------------------------------------------------------------------
    % CALCULATE ENERGY PARAMETERS
    %--------------------------------------------------------------------------
    
    Early5_80JLFC = permute(sum(abs(repmat(Early5_80(:,1,:),[1,3,1]) .* Early5_80(:,2:4,:))),[3,2,1]);
    Direct10 = permute(sum(Direct10.^2),[3,2,1]);
    Early5_80 = permute(sum(Early5_80.^2),[3,2,1]);
    Early0_80 = permute(sum(Early0_80.^2),[3,2,1]);
    Late80_end = permute(sum(Late80_end.^2),[3,2,1]);
    All_IR = permute(sum(All_IR.^2),[3,2,1]);
    
    
    JLF = pow2db(Early5_80(:,2:4)./repmat(Early0_80(:,1),[1,3]));
    
    JLFC = pow2db(Early5_80JLFC./repmat(Early0_80(:,1),[1,3]));
    
    out.funcallback.name = 'LateralFractions.m';
    out.funcallback.inarg = {fs, startthresh,threshmethod, bpo, omniweight, doplot};
    
    % OUTPUT TABLE
    out.tables = [];
    
    f = figure('Name','Spatial enegy fractions', ...
        'Position',[200 200 620 360]);
    %[left bottom width height]
    dat1 = [JLF(:,1)';JLFC(:,1)'];
    dat2 = [JLF(:,2)';JLFC(:,2)'];
    dat3 = [JLF(:,3)';JLFC(:,3)'];
    cnames1 = num2cell(bandfc);
    rnames1 = {'JLF (dB)', 'JLFC (dB)'};
    t1 =uitable('Data',dat1,'ColumnName',cnames1,'RowName',rnames1);
    t2 =uitable('Data',dat2,'ColumnName',cnames1,'RowName',rnames1);
    t3 =uitable('Data',dat3,'ColumnName',cnames1,'RowName',rnames1);
    set(t1,'ColumnWidth',{60});
    set(t2,'ColumnWidth',{60});
    set(t3,'ColumnWidth',{60});
    
    
    [~,tables] = disptables(f,[t1 t2 t3],{'Chan 2 - Energy Fraction';'Chan 3 - Energy Fraction';'Chan 4 - Energy Fraction'});
    out.tables = [out.tables tables];
    
    if doplot
        ymax = 10*ceil(max(max([JLF,JLFC]+5))/10);
        ymin = 10*floor(min(min([JLF,JLFC]+5))/10);
        figure('name', 'Spatial enegy fractions')
        for ch = 1:chans-1
            subplot(3,2,2*ch-1)
            width = 0.5;
            bar(1:length(bandfc),JLF(:,ch),width,'FaceColor',[1,0.3,0.3],...
                'EdgeColor',[0,0,0],'DisplayName', 'JLF','BaseValue',ymin);
            hold on
            
            % x-axis
            set(gca,'XTick',1:length(bandfc),'XTickLabel',num2cell(bandfc))
            if bpo==3
                xlabel('1/3-Octave Band Centre Frequency (Hz)')
            else
                xlabel('Octave Band Centre Frequency (Hz)')
            end
            
            % y-axis
            ylabel('JLF (dB)')
            ylim([ymin ymax])
            
            legend 'off'
            
            for k = 1:length(bandfc)
                text(k-0.25,JLF(k,ch)+2.5, ...
                    num2str(round(JLF(k,ch)*10)/10),'Color',[1,0.3,0.3])
            end
            
            
            
            subplot(3,2,2*ch)
            width = 0.5;
            bar(1:length(bandfc),JLFC(:,ch),width,'FaceColor',[0.2,0.6,0.2],...
                'EdgeColor',[0,0,0],'DisplayName', 'JLFC','BaseValue',ymin);
            hold on
            
            % x-axis
            set(gca,'XTick',1:length(bandfc),'XTickLabel',num2cell(bandfc))
            if bpo==3
                xlabel('1/3-Octave Band Centre Frequency (Hz)')
            else
                xlabel('Octave Band Centre Frequency (Hz)')
            end
            
            % y-axis
            ylabel('JLFC (dB)')
            ylim([ymin ymax])
            
            legend 'off'
            
            for k = 1:length(bandfc)
                text(k-0.25,JLFC(k,ch)+2.5, ...
                    num2str(round(JLFC(k,ch)*10)/10),'Color',[0.2,0.6,0.2])
            end
        end
        
    end
else
    out = [];
end
% eof

