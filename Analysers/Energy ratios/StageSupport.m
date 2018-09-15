function out = StageSupport(data, fs, startthresh,threshmethod, bpo, doplot)
% This function calculates stage support parameters from a room impulse 
% response (data).
%
% The start of the impulse response is detected using a decibel value
% specified by startthresh, and -20 dB is recommended by ISO3382-1.
%
% Normally start threshold detection and truncation should be done prior to
% filtering. The filterbanks used include zero-padding, to avoid loosing
% energy (from the filters' time response). However, if the impulse
% response has a highly variable group delay, it may be helpful to try
% start threshold detection and truncation after filtering. Actually it is
% better to avoid this problem in the first place by proper selection of
% measurement equipment.
%
% Calculations would normally be done in octave bands (bpo = 1), but
% 1/3-octave band analysis is also available.
%
% Currently this function only analyses the first two channels of
% multichannel audio input.
%
% Output values include ST Early (which is the same as ST1), ST2, and ST
% Late.
%
% Code by Densil Cabrera and Grant Cuthbert
% Version 1.01 (7 September 2014)


if nargin < 6, doplot = 1; end
if nargin < 5, bpo = 1; end
if nargin < 4, threshmethod = 0; end
if nargin < 3
    startthresh = -20;
    %dialog box for settings
    prompt = {'Threshold for IR start detection', ...
        'Detect start before filtering (0) or after filtering (1)',...
        'Bands per octave (1 | 3)', ...
        'Plot (0|1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {num2str(startthresh),'0','1','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        startthresh = str2num(answer{1,1});
        threshmethod = str2num(answer{2,1});
        bpo = str2num(answer{3,1});
        doplot = str2num(answer{4,1});
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
        warndlg('The StageSupport function is not suitable for multiband audio. Multiband audio has been mixed down, but results should be interpreted with caution.')
        IR = mean(IR,3);
    end

    if chans > 2
        IR = IR(:,1:2);
        warndlg('Only the first two channels are analysed')
        chans = 2;
    end

    if len/fs < 1
        warndlg('The impulse response is not sufficiently long for STLate. Zero padding has been applied, but results should be interpreted with caution.')
        IR = [IR; zeros(round(fs),chans)];
    end
    
    if threshmethod == 0
        % auto-crop start of individual audio columns
        IR = autocropstart_aarae(IR,startthresh,1);
        
        % Truncate sections of the IR prior to filtering
        Direct10 = IR(1:1+floor(fs*0.01),:); % Truncate first 10 ms
        Early20_100 = IR(1+floor(fs*0.02):1+floor(fs*0.1),:); % 20-100 ms
        Late100_200 = IR(1+floor(fs*0.1):1+floor(fs*0.2),:); % 100-200 ms
        All20_1000 = IR(1+floor(fs*0.02):1+floor(fs*1),:); % 20-1000 ms
        Late100_1000 = IR(1+floor(fs*0.1):1+floor(fs*1),:); % 100-1000 ms
    

    %--------------------------------------------------------------------------
    % FILTERING
    %--------------------------------------------------------------------------

    % Use inbuilt AARAE filterbank
    if bpo == 3
        bandfc = [100,125,160,200,250,315,400,500,630,800,1000,1250,1600, ...
            2000,2500,3150,4000,5000];
        All_IR = thirdoctbandfilter_viaFFT(IR,fs,bandfc,15,round(0.1*fs)); % IR
        Direct10 = thirdoctbandfilter_viaFFT(Direct10,fs,bandfc,15,round(0.1*fs)); 
        Early20_100 = thirdoctbandfilter_viaFFT(Early20_100,fs,bandfc,15,round(0.1*fs));
        Late100_200 = thirdoctbandfilter_viaFFT(Late100_200,fs,bandfc,15,round(0.1*fs));
        All20_1000 = thirdoctbandfilter_viaFFT(All20_1000,fs,bandfc,15,round(0.1*fs));
        Late100_1000 = thirdoctbandfilter_viaFFT(Late100_1000,fs,bandfc,15,round(0.1*fs));

    else
        bandfc = [125,250,500,1000,2000,4000];
        All_IR = octbandfilter_viaFFT(IR,fs,bandfc,12,round(0.1*fs)); % IR
        Direct10 = octbandfilter_viaFFT(Direct10,fs,bandfc,12,round(0.1*fs)); 
        Early20_100 = octbandfilter_viaFFT(Early20_100,fs,bandfc,12,round(0.1*fs));
        Late100_200 = octbandfilter_viaFFT(Late100_200,fs,bandfc,12,round(0.1*fs));
        All20_1000 = octbandfilter_viaFFT(All20_1000,fs,bandfc,12,round(0.1*fs));
        Late100_1000 = octbandfilter_viaFFT(Late100_1000,fs,bandfc,12,round(0.1*fs));

    end
    end
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    % TRUNCATION AFTER FILTERING (NOT RECOMMENDED, BUT MIGHT HELP IF GROUP
    % DELAY OF THE IR START VARIES A LOT)
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
        Early20_100 = IR(1+floor(fs*0.02):1+floor(fs*0.1),:,:); % 20-100 ms
        Late100_200 = IR(1+floor(fs*0.1):1+floor(fs*0.2),:,:); % 100-200 ms
        All20_1000 = IR(1+floor(fs*0.02):1+floor(fs*1),:,:); % 20-1000 ms
        Late100_1000 = IR(1+floor(fs*0.1):1+floor(fs*1),:,:); % 100-1000 ms
    end
    
    
    
    
    
    %--------------------------------------------------------------------------
    % CALCULATE ENERGY PARAMETERS
    %--------------------------------------------------------------------------


    Direct10 = permute(sum(Direct10.^2),[3,2,1]);
    Early20_100 = permute(sum(Early20_100.^2),[3,2,1]);
    Late100_200 = permute(sum(Late100_200.^2),[3,2,1]);
    Late100_1000 = permute(sum(Late100_1000.^2),[3,2,1]);
    All20_1000 = permute(sum(All20_1000.^2),[3,2,1]);

    ST1 = 10*log10(Early20_100 ./ Direct10); 
    ST2 = 10*log10(Late100_200 ./ Direct10); 
    STLate = 10*log10(Late100_1000 ./ Direct10);
    STAll = 10*log10(All20_1000 ./ Direct10);

    meanbands = find(bandfc>=200 & bandfc<=2500);
    ST1av = mean(ST1(meanbands,:),1);
    ST2av = mean(ST2(meanbands,:),1);
    STLateav = mean(STLate(meanbands,:),1);
    STAllav = mean(STAll(meanbands,:),1);


    %--------------------------------------------------------------------------
    % OUTPUT STRUCTURE (not used anymore)
    %--------------------------------------------------------------------------

%     if chans == 1
%         out.ST1 = [bandfc;ST1'];
%         out.ST2 = [bandfc;ST2'];
%         out.STLate = [bandfc;STLate'];
%         out.funcallback.name = 'StageSupport.m';
%         out.funcallback.inarg = {fs,startthresh,bpo};
% 
%     else %chans == 2
%         out.ST1_ch1 = [bandfc;ST1(:,1)'];
%         out.ST1_ch2 = [bandfc;ST1(:,2)'];
% 
%         out.ST2_ch1 = [bandfc;ST2(:,1)'];
%         out.ST2_ch2 = [bandfc;ST2(:,2)'];
% 
%         out.STLate_ch1 = [bandfc;STLate(:,1)'];
%         out.STLate_ch2 = [bandfc;STLate(:,2)'];
% 
%         out.funcallback.name = 'StageSupport.m';
%         out.funcallback.inarg = {fs,startthresh,bpo};
%     end % if chans

    out.funcallback.name = 'StageSupport.m';
    out.funcallback.inarg = {fs,startthresh,threshmethod,bpo};

    % OUTPUT TABLE
    
    
    out.tables = [];
    for ch = 1:chans
        f = figure('Name',['Stage Support, channel: ', num2str(ch)], ...
            'Position',[200 200 620 360]);
        %[left bottom width height]
        dat1 = [ST1(:,ch)';ST2(:,ch)';STLate(:,ch)';STAll(:,ch)'];
        cnames1 = num2cell(bandfc);
        rnames1 = {'ST1 or ST Early (dB)', 'ST2 (dB)', 'ST Late (dB)', 'ST All (dB)'};
        t1 =uitable('Data',dat1,'ColumnName',cnames1,'RowName',rnames1);
        set(t1,'ColumnWidth',{60});
        
        dat2 = [ST1av(:,ch)';ST2av(:,ch)';STLateav(:,ch)';STAllav(:,ch)'];
        cnames2 = {'Spectral Average (250 Hz - 2kHz Octave Bands)'};
        rnames2 = {'ST1 or ST Early (dB)','ST2 (dB)', ...
            'ST Late (dB)', 'ST All (dB)'};
        t2 =uitable('Data',dat2,'ColumnName',cnames2,'RowName',rnames2);
        
        [~,tables] = disptables(f,[t2 t1],{['Chan' num2str(ch) ' - Spectral average'],['Chan' num2str(ch) ' - Stage support']});
        out.tables = [out.tables tables];
    end
    if doplot
            ymax = 10*ceil(max(max([ST1,ST2,STLate]+5))/10);
            ymin = 10*floor(min(min([ST1,ST2,STLate]-5))/10);
        for ch = 1:chans
            figure('name', ['Stage Support, Channel ', ...
                num2str(ch)])
            
            
            subplot(3,1,1)
            width = 0.5;
            bar(1:length(bandfc),ST1(:,ch),width,'FaceColor',[1,0.3,0.3],...
                'EdgeColor',[0,0,0],'DisplayName', 'STEarly','BaseValue',ymin);
            hold on

            % x-axis
            set(gca,'XTick',1:length(bandfc),'XTickLabel',num2cell(bandfc))
            if bpo==3
                xlabel('1/3-Octave Band Centre Frequency (Hz)')
            else
                xlabel('Octave Band Centre Frequency (Hz)')
            end

            % y-axis
            ylabel('ST Early (dB)')
            ylim([ymin ymax])

            legend 'off'

            for k = 1:length(bandfc)
                text(k-0.25,ST1(k,ch)+2.5, ...
                    num2str(round(ST1(k,ch)*10)/10),'Color',[1,0.3,0.3])
            end
            
            
            
            subplot(3,1,2)
            width = 0.5;
            bar(1:length(bandfc),ST2(:,ch),width,'FaceColor',[0.2,0.6,0.2],...
                'EdgeColor',[0,0,0],'DisplayName', 'ST2','BaseValue',ymin);
            hold on

            % x-axis
            set(gca,'XTick',1:length(bandfc),'XTickLabel',num2cell(bandfc))
            if bpo==3
                xlabel('1/3-Octave Band Centre Frequency (Hz)')
            else
                xlabel('Octave Band Centre Frequency (Hz)')
            end

            % y-axis
            ylabel('ST2 (dB)')
            ylim([ymin ymax])

            legend 'off'

            for k = 1:length(bandfc)
                text(k-0.25,ST2(k,ch)+2.5, ...
                    num2str(round(ST2(k,ch)*10)/10),'Color',[0.2,0.6,0.2])
            end
            
            
            
            subplot(3,1,3)
            width = 0.5;
            bar(1:length(bandfc),STLate(:,ch),width,'FaceColor',[0.3,0.3,1],...
                'EdgeColor',[0,0,0],'DisplayName', 'STLate','BaseValue',ymin);
            hold on

            % x-axis
            set(gca,'XTick',1:length(bandfc),'XTickLabel',num2cell(bandfc))
            if bpo==3
                xlabel('1/3-Octave Band Centre Frequency (Hz)')
            else
                xlabel('Octave Band Centre Frequency (Hz)')
            end

            % y-axis
            ylabel('ST Late (dB)')
            ylim([ymin ymax])

            legend 'off'

            for k = 1:length(bandfc)
                text(k-0.25,STLate(k,ch)+2.5, ...
                    num2str(round(STLate(k,ch)*10)/10),'Color',[0.3,0.3,1])
            end
            
            
        end
        
    end
else
    out = [];
end
% eof

