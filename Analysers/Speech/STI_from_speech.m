function [Verbose,varargout] = STI_from_speech(in, fs, cal, ReferenceChannel, AuditoryMasking, doplot)
% This function estimates the speech transmission index from a 2-channel
% speech recording. One of the channels is used as reference (i.e., it is
% 'dry' speech without noise, which is input to the system being tested),
% while the other channel is the output of the system (e.g., with
% reverberation and noise).
%
% The modulation depth ratio of each channel is derived, from which the
% modulation transfer function is estimated. Following this, speech
% transmission index is calculated in the normal manner.
%
% One way of using this is to play a 'dry' recording of speech into the
% system, and record the output. Another way is to use real speech, with
% two microphones (one near the talker's mouth, and one at the measurement
% position).
%
% This code is only loosely based on a standard - it is mainly for
% exploration.
%
% Code by Densil Cabrera
% version 0 (22 October 2013) NOT TESTED PROPERLY!!


% INPUTS AND SETTINGS
if isstruct(in)
    audio = in.audio;
    fs = in.fs;
    if isfield(in,'cal')
        cal = in.cal;
        calgain = 10^(cal/20);
    else
        cal = 0;
        calgain = 1;
    end
else
    audio = in;
    if nargin < 3
        cal = inputdlg({'Calibration offset [dB]'},...
                           'Cal',1,{'0'});
        cal = str2num(char(cal));
    end
    calgain = 10^(cal/20);
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end
if nargin < 6, doplot = 0; end
if nargin < 5, AuditoryMasking = 1; end
if nargin < 4
    ReferenceChannel =1;
    % dialog box to get auditory masking, noise correction  and
    % reference channel settings
    prompt = {'Reference channel (1 | 2):', ...
        'Auditory Masking (0 | 1 | 2):', ...
        'Calibration offset (dB):', ...
        'Plotting (0 | 1):'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'1','1',num2str(round(10*20*log10(calgain))/10),'1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(answer)
        Verbose = [];
        return
    else
        ReferenceChannel = str2num(answer{1,1});
        AuditoryMasking = str2num(answer{2,1});
        calgain = 10.^(str2num(answer{3,1})/20);
        doplot = str2num(answer{4,1});
    end
end

if ~isempty(audio) && ~isempty(fs) && ~isempty(calgain) && ~isempty(ReferenceChannel) && ~isempty(AuditoryMasking) && ~isempty(doplot)
    audio = squeeze(mean(audio,3)); % mixdown 3rd dimension if it exists
    [len,chans] = size(audio);
    if chans < 2
        Verbose = [];
        warndlg('Two channels are required for STI_from_speech','AARAE info')
        return
    end
    audio = audio(:,1:2); % delete additional channels if they exist
    chans = 2;

    if ReferenceChannel == 2
        audio = circshift(audio,2); % swap channels
    end

    audio = audio .* calgain; % apply calibration


    % Nyquist frequency
    Nyquist=fs/2;

    % Time in seconds for each sample
    t=((1:len)-1)'./fs;


    % Define the octave band filter parameters for the basic filterbank
    bandnumber=21:3:39; % filter band numbers
    fc=10.^(bandnumber./10); % filter centre frequencies in Hz
    bandwidth = 1; % set to 0.5 instead of 1 if you want to try half-octave bandwidths
    f_low=fc./10^(0.15*bandwidth); % low cut-off frequency in Hz
    f_hi=fc.*10^(0.15*bandwidth); % high cut-off frequency in Hz

    % use filter order of 6 (half-order = 3)
    halforder = 3; % half of the filter order: a value of 3 yields a 6th order filter

    b = zeros(halforder*2+1,length(fc)); % pre-allocate filter coefficients
    a = b; % pre-allocate filter coefficients

    % calculate filter coefficients
    for k = 1:length(fc)
        [b(:,k), a(:,k)]=butter(halforder, [f_low(k)/Nyquist f_hi(k)/Nyquist]);
    end

    % Filter and square the data
    Intensity = zeros(len,chans,7);
    for k = 1:7
        Intensity(:,:,k)=(filter(b(:,k),a(:,k), audio)).^2; % filter & square
    end

    I = squeeze(mean(Intensity(:,2,:))); % mean square intensity of the received channel
    Level = 10*log10(I); % intensity expressed in decibels

    % low-pass filter at 100 Hz to filter the envelope (as per standard)
    [b, a] = butter(6,100/Nyquist,'low');
    Intensity = filter(b, a, Intensity);





    % CALCULATE MODULATION DEPTH RATIO, MDR
    % AND MODULATION TRANSFER FUNCTION, MTF

    % Modulation frequencies in Hz
    Fm = repmat(10.^((-2:11)/10),[7,1]);

    % number of whole number cycles to use for each modulation frequency
    Fm_cycles = floor(len .* Fm./fs);

    % number of samples to use for each modulation frequency
    Fm_len = floor(fs.*Fm_cycles ./ Fm);

    % time replicated across channels
    tch = repmat(t,1,chans);

    mdr = zeros(1,chans,7,14);
    for k = 1:7
        for j = 1:14
            % modulation depth of the reference and received signals
            mdr(1,:,k,j) = 2 * ((sum(Intensity(1:Fm_len(k,j),:,k) ...
                .* sin(2*pi*Fm(k,j)*tch(1:Fm_len(k,j),:)))).^2 ...
                + (sum(Intensity(1:Fm_len(k,j),:,k) ...
                .* cos(2*pi*Fm(k,j)*tch(1:Fm_len(k,j),:)))).^2).^0.5 ...
                ./ sum(Intensity(1:Fm_len(k,j),:,k));
        end

    end

    % Transform mdr to mtf using channel 1 as the reference
    MTF = mdr(1,2,:,:) ./ mdr(1,1,:,:);
    MTF = squeeze(MTF); % reduce to 2 dimensions




    if AuditoryMasking
        % AUDITORY MASKING

        amf=zeros(1,7);     % Auditory Masking Factor
        amfdB=zeros(1,7); % Auditory Masking Factor in dB
        Iam=zeros(1,7);   % Masking Intensity


        % auditory masking from 2011 version of the standard
        for k=2:7
            if Level(k-1)< 63
                amfdB(k)=0.5*Level(k-1)-65;
            elseif Level(k-1) < 67;
                amfdB(k)=1.8*Level(k-1)-146.9;
            elseif Level(k-1) < 100;
                amfdB(k)=0.5*Level(k-1)-59.8;
            else
                amfdB(k)=-10;
            end
            amf(k)=10^(amfdB(k)/10);
            Iam(k)=I(k-1)*amf(k);
        end

        % Absolute Speech Reception Threshold
        ART=[46 27 12 6.5 7.5 8 12];

        % Intensity of threshold
        Irt=10.^(ART./10);



        for k = 1:7
            MTF(k,:) = MTF(k,:) .* (repmat(I(k),[1,14])) ...
                ./ (repmat(I(k),[1,14]) ...
                + repmat(Iam(k),[1,14]) ...
                + Irt(k));
        end

    end



    MTF(MTF>1) = 1;


    % Effecive signal to noise ratio (SNR)
    SNReff=10*log10(MTF./(1-MTF));

    % limit values to -15 <= effSNR <= 15 dB
    SNReff(SNReff>15)=15;
    SNReff(SNReff<-15)=-15;

    % Calculate Transmission Index (TI)
    % and averaged Modulation Tranasfer Index (MTI)
    TI=((SNReff+15)./30)';
    MTI=mean(TI,1);


    % STI Male
    alpha=[0.085 0.127 0.230 0.233 0.309 0.224 0.173];
    beta=[0.085 0.078 0.065 0.011 0.047 0.095];
    MTI_alpha=alpha.*MTI;
    MTI_beta=beta.*sqrt(MTI(1:length(beta)).*MTI(2:length(beta)+1));
    M_STI=sum(MTI_alpha)-sum(MTI_beta);


    % STI Female
    alpha=[0 0.117 0.223 0.216 0.328 0.250 0.194];
    beta=[0 0.099 0.066 0.062 0.025 0.076];
    MTI_alpha=alpha.*MTI;
    MTI_beta=beta.*sqrt(MTI(1:length(beta)).*MTI(2:length(beta)+1));
    F_STI=sum(MTI_alpha)-sum(MTI_beta);


    % STIPA(IR) (2011 version)
    alpha=[0.085 0.127 0.230 0.233 0.309 0.224 0.173];
    beta=[0.085 0.078 0.065 0.011 0.047 0.095];
    MTI_STIPA = [mean([TI(5,1); TI(12,1)]), ...
        mean([TI(3,2); TI(10,2)]), ...
        mean([TI(1,3); TI(8,3)]), ...
        mean([TI(6,4); TI(13,4)]), ...
        mean([TI(4,5); TI(11,5)]), ...
        mean([TI(2,6); TI(9,6)]), ...
        mean([TI(7,7); TI(14,7)])];
    MTI_alpha=alpha.*MTI_STIPA;
    MTI_beta=beta.*sqrt(MTI_STIPA(1:length(beta)) ...
        .*MTI_STIPA(2:length(beta)+1));
    STIPA=sum(MTI_alpha)-sum(MTI_beta);


    MTF=MTF';

    Verbose.M_STI = M_STI;
    Verbose.F_STI = F_STI;
    Verbose.MTF = MTF;
    Verbose.MTI = MTI;
    Verbose.Level = Level;
    Verbose.funcallback.name = 'STI_from_speech.m';
    Verbose.funcallback.inarg = {fs,cal,ReferenceChannel,AuditoryMasking,doplot};
    varargout{1} = M_STI;
    varargout{2} = F_STI;

    %***************************************************************
    % Plotting

    if doplot


        % Create figure
        figure('Name', 'STI estimated from speech');

        % plot of mtfs
        subplot(3,1,1)
        hold on

        % x axis
        set(gca,'XTickLabel',{'0.63', '0.8', '1', '1.25', '1.6', '2',...
            '2.5', '3.15','4','5','6.3','8','10','12.5'}, ...
            'XTick',1:14)
        xlim([1 14]);

        % rainbow colours
        colors = [255, 0, 0; ... % red
            255, 128, 0; ... % orange
            204, 204, 0; ... % dark yellow
            0, 204, 0; ... % mid green
            0, 204, 204; ... % dark cyan
            0, 0, 255; ... % blue
            127, 0, 255]; % violet
        colors = colors / 255; % rescale to 0-1 range

        % plot the mtf for each octave band (DisplayName is used if a
        % legend is added). Use rainbow colours as defined above.

        plot(MTF(:,7),'Color',colors(7,:),'DisplayName','8 kHz');
        plot(MTF(:,6),'Color',colors(6,:),'DisplayName','4 kHz');
        plot(MTF(:,5),'Color',colors(5,:),'DisplayName','2 kHz');
        plot(MTF(:,4),'Color',colors(4,:),'DisplayName','1 kHz');
        plot(MTF(:,3),'Color',colors(3,:),'DisplayName','500 Hz');
        plot(MTF(:,2),'Color',colors(2,:),'DisplayName','250 Hz');
        plot(MTF(:,1),'Color',colors(1,:),'DisplayName','125 Hz');

        % legend
        legend('show','Location','EastOutside');

        % y-axis limits
        ylim([0 1]);

        % Create xlabel
        xlabel('Modulation Frequency (Hz)');

        % Create ylabel
        ylabel('Modulation Transfer Coefficient');

        % Create title

            title(['Male STI ', num2str(round(M_STI*100)/100), ...
                '; Female STI ', num2str(round(F_STI*100)/100), ...
                '  (AM=', num2str(AuditoryMasking),')'])

        hold off

        subplot(3,1,2)
        % Bar plot of modulation transfer indices

        % Plot the values
        bh=bar(MTI);

        % The following code enables each bar to be coloured individually,
        % using the same colours as the mtf subplot (hence the bar plot
        % acts as a legend for the mtf plot).
        chil=get(bh,'children');
        cd=repmat(1:7,5,1);
        cd=[cd(:);nan];
        colormap(colors);
        set(chil,'facevertexcdata',cd);

        % x-axis
        set(gca,'XTickLabel',{'125', '250', '500', '1k', '2k', '4k', '8k'})
        xlabel('Octave Band Centre Frequency (Hz)')

        % y-axis
        ylabel('Modulation Transfer Index')
        ylim([0 1])

        % show the value of each bar
        for k = 1:7
            if MTI(k) < 0.8
                % black number above the bar
                text(k-0.25,MTI(k)+0.05, ...
                    num2str(round(MTI(k)*1000)/1000),'Color','k')
            else
                % or white number below the top of the bar
                % (for high MTI values)
                text(k-0.25,MTI(k)-0.05, ...
                    num2str(round(MTI(k)*1000)/1000),'Color',[1 1 1])
            end
        end

        subplot(3,1,3)
        % Plot of octave band signal and noise levels
        hold on

        % x-axis
        set(gca,'XTickLabel',{'125', '250', '500', '1k', '2k', '4k', '8k'})
        xlabel('Octave Band Centre Frequency (Hz)')

        % y-axis
        if max(Level)>30
        ymax = 10*ceil(max(Level)/10);
        ylim([0 ymax])
        end
        ylabel('SPL (dB)')

        % Sound pressure level of the received signal (& noise)
        plot(Level,'b','Marker','o','DisplayName','Received')


        if AuditoryMasking > 0

            % Sound pressure level of masking from the band below
            plot(10*log10(Iam),'Color',[0 0.5 0], ...
                'Marker','x','DisplayName','Masking')

            % Sound pressure level of the auditory reception threshold
            plot(10*log10(Irt),'k','LineStyle',':', ...
                'DisplayName','Threshold')

            % Sum of all sources of noise
            Isum = I' + Irt + Iam;
            plot(10*log10(Isum), 'r', ...
                'LineWidth',1,'LineStyle','--', ...
                'Marker','+', 'DisplayName','Total S+N')
            
            if isstruct(in)
                doresultleaf([Level,10*log10(Iam'),10*log10(Irt'),10*log10(Isum')],'SPL [dB]',{'Frequency'},...
                             'Frequency',            num2cell([125,250,500,1000,2000,4000,8000]), 'Hz', true,...
                             'Level',     {'Received SPL','Masking','Threshold','Total S+N'}, 'categorical', [],...
                             'name','Band_SPL');
            end
        end
        legend('show','Location','EastOutside');
        hold off
    end % if doplot
    if isstruct(in)
        mf = [0.63,0.8,1,1.25,1.6,2,2.5,3.15,4,5,6.3,8,10,12.5];
        doresultleaf(MTF,'Coefficient',{'Modulation_frequency'},...
                     'Modulation_frequency', num2cell(mf),                                'Hz', true,...
                     'Frequency',            num2cell([125,250,500,1000,2000,4000,8000]), 'Hz', false,...
                     'name','Modulation_TF_coef');
        if AuditoryMasking > 0
            doresultleaf([Level,10*log10(Iam'),10*log10(Irt'),10*log10(Isum')],'SPL [dB]',{'Frequency'},...
                         'Frequency', num2cell([125,250,500,1000,2000,4000,8000]),        'Hz',          true,...
                         'Level',     {'Received SPL','Masking','Threshold','Total S+N'}, 'categorical', [],...
                         'name','Band_SPL');
        else
            doresultleaf([Level,10*log10(Iam'),10*log10(Irt'),10*log10(Isum')],'SPL [dB]',{'Frequency'},...
                         'Frequency', num2cell([125,250,500,1000,2000,4000,8000]), 'Hz',          true,...
                         'Level',     {'Received SPL'},                            'categorical', [],...
                         'name','Band_SPL');
        end
    end
end
% eof