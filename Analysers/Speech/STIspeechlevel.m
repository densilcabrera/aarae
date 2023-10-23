function [OUT,varargout] = STIspeechlevel(speech, fs, cal, doplot, wintime, caltone)
% This function calculates the operational speech level of a recording of
% speech as described by Annex J of IEC 60268-16 (2011), to support the
% calculation of speech transmission index (STI).
%
% Code by Densil Cabrera
% version 1.01 (16 November 2013)


if isstruct(speech)
    speech = choose_from_higher_dimensions(speech,2,1);
    signal = speech.audio;
    fs = speech.fs;
    if isfield(speech,'cal')
        cal = speech.cal;
    else
        cal = zeros(1,size(signal,2));
    end
    if isfield(speech,'chanID')
        chanID = speech.chanID;
    else
        chanID = cellstr([repmat('Chan',size(signal,2),1) num2str((1:size(signal,2))')]);
    end
else
    signal = speech;
    if nargin < 3
        cal = inputdlg({'Calibration offset [dB]'},...
            'Cal',1,{'0'});
        cal = str2num(char(cal));
    end
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
            'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end
%********************************************
% Calibration of the input signal

if nargin > 5 && ~isempty(caltone)
    % In this case, a calibration tone has been provided as input argument
    % 6
    % assume caltone = 1 kHz, and A-weight to filter out lf noise
    caltone = Aweight(caltone(:,1), fs);

    Amplitude = mean(caltone.^2).^0.5;
    Acal = 10^(cal/20);
    gain = Acal / Amplitude;
    if doplot
        figure
        plot(((1:length(caltone))-1)./fs, ...
            10*log10((caltone*gain).^2))
        hold on
        plot([0 (length(caltone)-1)/fs], [cal cal], ...
            'linestyle', '--', 'color', [0.5, 0.5, 0.5]);
        xlabel('Time (s)')
        title(['Calibration tone, LAeq ', num2str(cal), ' dBA'])
        ylabel('Squared calibration tone, in dB')
        hold off
    end
elseif  ~exist('cal','var')
    % default calibration factor
    caltone = [];
    cal = 0;
    gain = 1;
else
    % user-defined arbitrary calibration, using input argument 3
    caltone = [];
    gain = 10.^(cal./20);
end
% default window length, in milliseconds
if nargin < 5, wintime = 15; end
if nargin < 4, doplot = 1; end
if nargin < 3
    prompt = {'Calibration offset (dB):', ...
        'Window time (ms):', ...
        'Make Plot (0 | 1):'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {num2str(20*log10(gain)),'15','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if isempty(answer)
        OUT = [];
        return
    else
        cal = str2num(answer{1,1});
        gain = 10.^(cal./20);
        wintime = str2num(answer{2,1});
        doplot = str2num(answer{3,1});
    end
end

if ~isempty(signal) && ~isempty(fs) && ~isempty(cal) && ~isempty(doplot) && ~isempty(wintime)
    % remove 3rd dimension if it exists
    signal = squeeze(signal(:,:,1));
    [len,chans] = size(signal);

    % calibration adjustment of signal
    if length(gain) == chans
        signal = signal .* repmat(gain,[len 1]);
    else
        signal = singal .* gain;
    end

    % window length in samples
    winlength = round(wintime*fs/1000);

    % window overlap
    overlap = 0.9;
    step = round(winlength*(1-overlap));
    nwin =floor((len-winlength)/step);


    if doplot
        disp(['Number of windows: ', num2str(nwin)])
    end

    % ************************************
    % A-weight the speech recording
    signal1 = Aweight(signal,fs);

    % determine the sound intensity in each window
    I = zeros(nwin,chans);
    for n = 1:nwin
        a = (n - 1) * step + 1 ;
        b = a + winlength - 1;
        I(n,:) = mean(signal1(a:b,:).^2);
    end

    LL.original = 10*log10(I);

    Isort = sort(I,1,'descend');

    LL.sort = 10*log10(Isort);
    % cumulative mean, in decibels
    LL.cummean = 10*log10(cumsum(Isort) ./ repmat((1:length(Isort))',[1,chans]));

    ind = zeros(1,chans);
    L.speechJ2 = zeros(1,chans);

    for ch = 1:chans
        try
            ind(ch) = find(LL.sort(:,ch) <= LL.cummean(:,ch)-14,1,'first');
            L.speechJ2(ch) = LL.cummean(ind(ch),ch);
        catch
            L.speechJ2(ch) = nan;
        end
    end

    % Equivalent A-weignted sound pressure level
    L.eq = 10*log10(mean(I));

    % Simple estimate of speech level, by adding 3 dB to LAeq, based on
    % IEC 60268-16 (2011) J.4
    L.speechJ4 = L.eq + 3;

    % Maximum A-weighted sound pressure level
    L.max = LL.sort(1,:);

    % 95th percentile of A-weighted sound pressure level, which by convention
    % is called L5 in acoustics. Note that this value is calculated using the
    % windows specified in this function, rather than by using conventional
    % 'fast' integration
    L.p5 = prctile(LL.sort, 95);

    % 90th percentile, L10
    L.p10 = prctile(LL.sort, 90);

    % Median A-weighted sound pressure level (L50)
    L.median = median(LL.sort);

    % 10th percentile, L90
    L.p90 = prctile(LL.sort,10);

    % 5th percentile, L5
    L.p95 = prctile(LL.sort,5);

    % minimum A-weighted sound pressure level
    L.min = LL.sort(end,:);


    %*********************************************




    signal2 = zeros(length(signal),chans);
    for ch = 1:chans
        onoff = I(:,ch) >= I(ind(ch));
        winind = round(((1:nwin))'*step -(step-1)).* onoff;


        for n = 1:nwin
            if winind(n) > 0
                signal2(winind(n):winind(n)+winlength-1, ch) = ...
                    signal(winind(n):winind(n)+winlength-1, ch);
            end
        end

    end
    fc = [125,250,500,1000,2000,4000,8000];
    filterorder = 24; % very steep filter skirts improve the performance
    P_octave = octbandfilter_viaFFT(signal2,fs,...
        fc,filterorder);

    % OLD OCTAVE BAND FILTERS - DO NOT USE
    % % Define the octave band filter parameters
    % Nyquist = fs/2;
    % bandnumber=21:3:39; % filter band numbers
    % fc=10.^(bandnumber./10); % filter centre frequencies in Hz
    % bandwidth = 1; % set to 0.5 instead of 1 if you want to try half-octave bandwidths
    % f_low=fc./10^(0.15*bandwidth); % low cut-off frequency in Hz
    % f_hi=fc.*10^(0.15*bandwidth); % high cut-off frequency in Hz
    % halforder = 3; % half of the filter order: a value of 3 yields a 6th order filter
    %
    % P_octave = zeros(length(signal2), chans, length(fc));
    % for k=1:length(fc);
    %     % use filter order of 6 (half-order = 3)
    %     [b, a]=butter(halforder, [f_low(k)/Nyquist f_hi(k)/Nyquist]);
    %     P_octave(:,:,k)=filtfilt(b,a, signal2); % linear phase filter
    %
    % end
    onoffratio = zeros(1,chans);
    for ch = 1:chans
        onoffratio(1,ch) = ind(ch) / length(LL.sort);
    end
    % Define the octave band filter parameters
    % Nyquist = fs/2;
    % bandnumber=21:3:39; % filter band numbers
    % fc=10.^(bandnumber./10); % filter centre frequencies in Hz
    % bandwidth = 1; % set to 0.5 instead of 1 if you want to try half-octave bandwidths
    % f_low=fc./10^(0.15*bandwidth); % low cut-off frequency in Hz
    % f_hi=fc.*10^(0.15*bandwidth); % high cut-off frequency in Hz
    % halforder = 3; % half of the filter order: a value of 3 yields a 6th order filter
    % 
    % P_octave = zeros(length(signal2), chans, length(fc));
    % for k=1:length(fc);
    %     % use filter order of 6 (half-order = 3)
    %     [b, a]=butter(halforder, [f_low(k)/Nyquist f_hi(k)/Nyquist]);
    %     P_octave(:,:,k)=filtfilt(b,a, signal2); % linear phase filter
    % 
    % end
    % onoffratio = zeros(1,chans);
    % for ch = 1:chans
    %     onoffratio(1,ch) = ind(ch) / length(LL.sort);
    % end



    Lspeechoct = 10 * log10(mean(P_octave.^2) ...
        ./ repmat(onoffratio,[1,1,length(fc)]));

    % A-weight and sum the octave band values to find an approximate A-weighted
    % speech level. This is done as a check, because the value should be
    % approximately the same as L.speechJ2 (and L.speechJ4)
    L.Aspeechsumoctbandscheck = 10*log10( 10.^((Lspeechoct(1,:,1)-16.1)./10) ...
        + 10.^((Lspeechoct(1,:,2)-8.6)./10) ...
        + 10.^((Lspeechoct(1,:,3)-3.2)./10) ...
        + 10.^(Lspeechoct(1,:,4)./10) ...
        + 10.^((Lspeechoct(1,:,5)+1.2)./10) ...
        + 10.^((Lspeechoct(1,:,6)+1)./10) ...
        + 10.^((Lspeechoct(1,:,7)-1.1)./10));

    OUT = speech;
    OUT.audio = signal2;
    OUT.L = L;
    OUT.LL = LL;
    OUT.funcallback.name = 'STIspeechlevel.m';
    OUT.funcallback.inarg = {fs,cal,doplot,wintime,caltone};
    varargout{1} = L;
    varargout{2} = LL;

    if isstruct(speech)
        doresultleaf(10*log10(I),'SPL [dB]',{'Time'},...
            'Time',    num2cell(((1:nwin)-1) * step / fs), 'Hz',          true,...
            'Channel', chanID,                             'categorical', [],...
            'name','Speech_level');
    end

    if doplot
        times = ((1:nwin)-1) * step / fs;
        ymin = floor(min(LL.sort(end,:))/10) *10;
        ymax = ceil(max(LL.sort(1,:))/10)*10;
        for ch = 1:chans

            figure('Name',['Operational Speech Level, Channel ', num2str(ch)])

            subplot(3,1,1)

            plot(times, 10*log10(I(:,ch)));
            title(['Channel ', num2str(ch)])
            xlabel('Time (s)')
            ylabel('SPL (dB)')
            ylim([ymin ymax])
            hold on
            plot([times(1) times(end)], [LL.sort(ind(ch),ch) LL.sort(ind(ch),ch)], ...
                'color', [0.5 0.5 0.5])
            plot([times(1) times(end)], [L.speechJ2(ch) L.speechJ2(ch)], ...
                'LineStyle','--', 'color', [0.5 0.5 0.5])
            hold off


            subplot(3,1,2)

            plot(LL.sort,'k')
            ylim([ymin ymax])
            title(['Speech level (J2) ', num2str(round(L.speechJ2(ch)*10)/10), ' dB(A)'])
            hold on
            plot(LL.cummean, 'r')
            plot([ind(ch) ind(ch)], [ymin ymax], 'color', [0.5 0.5 0.5])


            xlabel('Number of windows')
            ylabel('SPL (dB)')


            hold off
            subplot(3,1,3)
            bar(squeeze(Lspeechoct(1,ch,:)))
            set(gca,'XTickLabel',{'125', '250', '500', '1k', '2k', '4k', '8k'})
            xlabel('Octave Band Centre Frequency (Hz)')
            ylabel('Speech SPL (dB)')

            for k = 1:7
                % numbers on chart
                text(k-0.2,Lspeechoct(1,ch,k)-5, ...
                    num2str(round(Lspeechoct(1,ch,k)*10)/10),'Color','y')
            end

        end

        % Loop for replaying, saving and finishing
        %        choice = 0;
        %        wave = signal2./max(max(abs(signal2)));
        % loop until the user presses the 'Done' button
        %        while choice < 5
        %            choice = menu('What next?', ...
        %                'Table of Results', ...
        %                'Play original', ...
        %                'Play thresholded', ...
        %                'Save wav file', ...
        %                'Done');
        %            switch choice

        %                case 1
        f = figure('Name','Operational Speech Level', ...
            'Position',[200 200 440+60*chans 360]);
        %[left bottom width height]
        Loctband = permute(Lspeechoct(1,:,:),[3,2,1]);
        dat = [L.speechJ2;L.eq;L.speechJ4;L.max;L.p5;L.p10;L.median;L.p90; ...
            L.p95;L.min;Loctband];
        cnames = {'SPL (dB)'};
        rnames = {'A-Weighted Operational Speech Level J2', ...
            'Equivalent A-weighted SPL, LA,eq', ...
            'A-Weighted Operational Speech Level J4', ...
            'Maximum A-weighted Level, LA,max', ...
            '95th Percentile A-weighted, SPL LA,5', ...
            '90th Percentile A-weighted, SPL LA,10', ...
            'Median A-weighted SPL, LA,50', ...
            '10th Percentile A-weighted SPL, LA,90', ...
            '5th Percentile A-weighted SPL, LA,95', ...
            'Minimum A-weighted SPL, LA,min', ...
            'Octave Band SPL J2 125 Hz', ...
            'Octave Band SPL J2 250 Hz', ...
            'Octave Band SPL J2 500 Hz', ...
            'Octave Band SPL J2 1 kHz', ...
            'Octave Band SPL J2 2 kHz', ...
            'Octave Band SPL J2 4 kHz', ...
            'Octave Band SPL J2 8 kHz'};
        t =uitable('Parent',f,'Data',dat,'ColumnName',cnames,...
            'RowName',rnames,'Position',[20 20 400+60*chans 340]);
        %set(t,'ColumnWidth',{60});
        disptables(f,t)

        %                case 2
        %                    wave = signal./max(max(abs(signal)));
        %                    sound(wave, fs)
        %                case 3
        %                    wave = signal2./max(max(abs(signal2)));
        %                    sound(wave, fs)

        %                case 4
        %                    [filename, pathname] = uiputfile({'*.wav'},'Save as');
        %                    if ischar(filename)
        %                        audiowrite([pathname,filename], wave, fs);
        %                    end
        %            end
        %        end

    end
else
    OUT = [];
end % eof
end

function y = Aweight(x, fs)
% Filters input x with A-weighting and returns output y.

% THE FOLLOWING IS NOT COMPATIBLE WITH CURRENT MATLAB VERSION
% MATLAB Code
% Generated by MATLAB(R) 8.1 and the DSP System Toolbox 8.4.
% Generated on: 24-Aug-2013 14:55:39

% persistent Hd;
%
% if isempty(Hd)
%
%     WT    = 'A';    % Weighting type
%     Class = 1;      % Class
%     if nargin < 2, fs    = 48000; end  % Sampling Frequency
%
%     h = fdesign.audioweighting('WT,Class', WT, Class, fs);
%
%      Hd = design(h, 'ansis142', ...
%          'SOSScaleNorm', 'Linf');
% end
AWeighting = weightingFilter('A-weighting',fs);
y = AWeighting(x);

end