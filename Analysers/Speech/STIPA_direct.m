function out = STIPA_direct(in, fs, cal, refsignal, AuditoryMasking, doplot, dorefsig)
% This function analyses a previously recorded STIPA signal, and derives
% the value of STIPA, as described by IEC60268-16 (2011)
%
% Before running this function, play and record the STIPA test signal
% through the acoustic or audio system that you wish to test.
% Note that both the loudspeaker (if used) and microphone (if used) will
% need to be calibrated to make a correct measurement. The loudspeaker may
% also need to be equalized to produce the correct STIPA spectrum levels.
% Analyse the recorded signal using this function. The recommended signal
% duration is 20 s or greater (accuracy should improve with greater
% duration).
%
% This function gives the option of using a reference signal - which could
% either be the original STIPA signal (to correct for minor errors in the
% calculated modulation depth ratio), or the STIPA signal re-recorded
% in low-noise low-reverberant energy conditions (e.g., at the mouth
% reference point of a head and torso simulator) to correct for transducer
% effects. 
%
% This version uses AARAE's linear (actually zero) phase octave band filters. 
% Steep filter skirts return the best results, so currently 144 dB/oct 
% skirts are used (achieved by setting filterorder to 24).
%
% Code by Densil Cabrera
% version 1.04 (31 January 2014)


% INPUTS AND SETTINGS
if isstruct(in)
    in = choose_from_higher_dimensions(in,2,1);
    audio = in.audio;
    fs = in.fs;
    if isfield(in,'cal')
        cal = in.cal;
        calgain = 10.^(cal/20);
    else
        cal = 0;
        calgain = 1;
    end
    if isfield(in,'name') % Get the AARAE name if it exists
        name = in.name; % this is a string
    else
        name = ''; % empty string - can be concatenated without any problems
    end
else
    audio = in;
    if nargin < 3
        cal = inputdlg({'Calibration offset [dB]'},...
                           'Cal',1,{'0'});
        cal = str2num(char(cal));
    end
    calgain = 10.^(cal/20);
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
    name = '';
end

if nargin < 7, dorefsig = 0; end
if nargin < 6, doplot = 1; end
if nargin < 5, AuditoryMasking = 1; end
if nargin < 4
    refsignal = [];
    % dialog box to get auditory masking, noise correction  and
    % auralization settings
    prompt = {'Auditory Masking (0 | 1):', ...
        'Calibration offset (dB):', ...
        'Plotting (0 | 1):',...
        'Additional reference signal (0 | 1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'1',num2str(round(10*20*log10(calgain))/10),'1','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(answer)
        out = [];
        return
    else
        AuditoryMasking = str2num(answer{1,1});
        calgain = 10.^(str2num(answer{2,1})/20);
        doplot = str2num(answer{3,1});
        dorefsig = str2num(answer{4,1});
    end
end

if ~isempty(audio) && ~isempty(fs) && ~isempty(cal) && ~isempty(AuditoryMasking) && ~isempty(doplot) && ~isempty(dorefsig)

    audio = mean(audio,3); % mixdown 3rd dimension if it exists
    chans = size(audio,2);
    audio = audio .* calgain; % apply calibration
    if isfield(in,'chanID')
        chanID = in.chanID;
    else
        chanID = cellstr([repmat('Chan',size(audio,2),1) num2str((1:size(audio,2))')]);
    end

    % Define the octave band filter parameters
    %bandnumber=21:3:39; % filter band numbers
    %fc=10.^(bandnumber./10); % filter centre frequencies in Hz
    fc = [125, 250, 500, 1000, 2000, 4000, 8000]; % nominal centre freq

    % list of modulation frequencies
    Fm = [1.6, 8;...
        1, 5;...
        0.63, 3.15;...
        2, 10;...
        1.25, 6.25;...
        0.8, 4;...
        2.5, 12.5];

    % -------------------------------------------------------------------------
    % MTF CALCULATIONS
    % -------------------------------------------------------------------------

    if dorefsig && isempty(refsignal)
        % Use a menu & dialog box to select a wav file or audio within AARAE
        if isstruct(in)
            selection = choose_audio; % call AARAE's choose_audio function
        end
        if ~isempty(selection) && isstruct(in)
            refsignal = selection.audio; % additional audio data
            fs2 = selection.fs; % sampling rate

            if ~(fs2 == fs)
                % match sampling rates if desired
                refsignal = resample(refsignal,fs,fs2);
            end
            [~, chans2, bands2] = size(refsignal); % new wave dimensions
            if bands2 > 1, refsignal = sum(refsignal,3); end % if multiband, mixdown
            if chans2>1, refsignal = refsignal(:,1); end % only 1 channel allowed for ref

            % CALCULATE MTF OF REFERENCE SIGNAL
            MTFref = STIPA_direct_MTF(refsignal,fs,fc,Fm);
        elseif ~isstruct(in)
            MTFref = STIPA_direct_MTF(refsignal,fs,fc,Fm);
        else
            dorefsig = 0;
        end
    elseif ~isempty(refsignal)
        MTFref = STIPA_direct_MTF(refsignal,fs,fc,Fm);
    end

    % CALCULATE MTF
    [MTF, I, Level] = STIPA_direct_MTF(audio,fs,fc,Fm);

    % adjust MTF by reference signal MTF
    if dorefsig
        MTF = MTF ./ MTFref;
    end

    % -------------------------------------------------------------------------
    % APPLY AUDITORY MASKING AND THRESHOLD
    % -------------------------------------------------------------------------

    if AuditoryMasking
        % AUDITORY MASKING
        % preallocate
        amf=zeros(1,chans,7);     % Auditory Masking Factor
        amfdB=zeros(1,chans,7); % Auditory Masking Factor in dB
        Iam=zeros(1,chans,7);   % Masking Intensity

        % auditory masking from 2011 version of the standard
        for ch = 1:chans
            for k=2:7
                if Level(1,ch,k-1)< 63
                    amfdB(1,ch,k)=0.5*Level(1,ch,k-1)-65;
                elseif Level(1,ch,k-1) < 67
                    amfdB(1,ch,k)=1.8*Level(1,ch,k-1)-146.9;
                elseif Level(1,ch,k-1) < 100
                    amfdB(1,ch,k)=0.5*Level(1,ch,k-1)-59.8;
                else
                    amfdB(1,ch,k)=-10;
                end
                amf(1,ch,k)=10^(amfdB(1,ch,k)/10);
                Iam(1,ch,k)=I(1,ch,k-1)*amf(1,ch,k);
            end
        end

        % Absolute Speech Reception Threshold
        ART=[46 27 12 6.5 7.5 8 12];

        % Intensity of threshold
        Irt=10.^(ART./10);

        for k = 1:7
            MTF(1,:,k,:) = MTF(1,:,k,:) .* (repmat(I(1,:,k),[1,1,1,2])) ...
                ./ (repmat(I(1,:,k),[1,1,1,2]) ...
                + repmat(Iam(1,:,k),[1,1,1,2]) ...
                + Irt(k));
        end
    end

    % -------------------------------------------------------------------------
    % FINAL STAGES OF STIPA CALCULATION
    % -------------------------------------------------------------------------

    % limiting the MTF values to 1 at the maximum.
    MTF(MTF>1) = 1;
    MTF = permute(MTF,[2,3,4,1]);

    % Effecive signal to noise ratio (SNR)
    SNReff=10*log10(MTF./(1-MTF));

    % limit values to -15 <= effSNR <= 15 dB
    SNReff(SNReff>15)=15;
    SNReff(SNReff<-15)=-15;

    % Calculate Transmission Index (TI)
    % and averaged Modulation Tranasfer Index (MTI)
    TI=(SNReff+15)./30;
    MTI=mean(TI,3);

    % STIPA (2011 version)
    alpha=repmat([0.085 0.127 0.230 0.233 0.309 0.224 0.173],[chans,1]);
    beta=repmat([0.085 0.078 0.065 0.011 0.047 0.095],[chans,1]);
    MTI_alpha=alpha.*MTI;
    MTI_beta=beta.*sqrt(MTI(:,1:6) ...
        .*MTI(:,2:7));
    STIPA=sum(MTI_alpha,2)-sum(MTI_beta,2);

    % -------------------------------------------------------------------------
    % CREATE OUTPUT STRUCTURE
    % -------------------------------------------------------------------------
    out.STIPA = STIPA;
    out.MTI = MTI;
    out.MTF = MTF;
    out.Level = Level;
    out.funcallback.name = 'STIPA_direct.m';
    out.funcallback.inarg = {fs,cal,refsignal,AuditoryMasking,doplot,dorefsig};

    
    % -------------------------------------------------------------------------
    % CREATE OUTPUT TABLES
    % -------------------------------------------------------------------------
    out.tables = [];
    for ch = 1:chans
        fig1 = figure('Name',['My results table for ' name]);
        table1 = uitable('Data',[permute(Level(:,ch,:),[2,3,1]); MTI(ch,:); permute(MTF(ch,:,:),[3,2,1])],...
            'ColumnName',{'125','250','500','1000','2000','4000','8000'},...
            'RowName',{'Level','MTI','MTF1','MTF2'});
        table2 = uitable('Data',STIPA(ch),...
            'ColumnName', {['Chan ' num2str(ch)]},...
            'RowName', {'STIPA'});
        [~,table] = disptables(fig1,[table1,table2]);
        out.tables = [out.tables table];
    end
    
    
    
    % -------------------------------------------------------------------------
    % PLOTTING
    % -------------------------------------------------------------------------

    if doplot
        for ch = 1:chans
            figure('Name', ['Channel: ', num2str(ch),' STIPA: ' num2str(STIPA(ch))])


            subplot1 = subplot(3,1,1);
            MTF_all = zeros(8,15);
            MTF_all(1,5) = MTF(ch,1,1);
            MTF_all(1,12) = MTF(ch,1,2);
            MTF_all(2,3) = MTF(ch,2,1);
            MTF_all(2,10) = MTF(ch,2,2);
            MTF_all(3,1) = MTF(ch,3,1);
            MTF_all(3,8) = MTF(ch,3,2);
            MTF_all(4,6) = MTF(ch,4,1);
            MTF_all(4,13) = MTF(ch,4,2);
            MTF_all(5,4) = MTF(ch,5,1);
            MTF_all(5,11) = MTF(ch,5,2);
            MTF_all(6,2) = MTF(ch,6,1);
            MTF_all(6,9) = MTF(ch,6,2);
            MTF_all(7,7) = MTF(ch,7,1);
            MTF_all(7,14) = MTF(ch,7,2);

            pcolor(MTF_all);
            colormap('hot')
            MTF_out(1:7,1:14,ch) = MTF_all(1:7,1:14);
            % x-axis
            set(gca,'XTickLabel',{'0.63', '0.8', '1', '1.25', '1.6', '2',...
                '2.5', '3.15','4','5','6.3','8','10','12.5'}, ...
                'XTick',1.5:1:14.5)
            xlim([1 15]);
            xlabel('Modulation Frequency (Hz)')

            % y-axis
            set(gca,'YTickLabel',{'125', '250', '500', '1k', '2k', '4k', '8k'},...
                'YTick',1.5:1:7.5)
            ylabel('Band (Hz)')
            ylim([1 8])

            colorbar('peer',subplot1);
            title('Sparse modulation transfer function')


            subplot(3,1,2)

            bar(squeeze(MTI(ch,:)),'r');

            % x-axis
            set(gca,'XTickLabel',{'125', '250', '500', '1k', '2k', '4k', '8k'})
            xlabel('Octave Band Centre Frequency (Hz)')

            % y-axis
            ylabel('MTI')
            ylim([0 1])

            % show the value of each bar
            for k = 1:7
                if MTI(ch,k) < 0.8
                    % black number above the bar
                    text(k-0.25,MTI(ch,k)+0.05, ...
                        num2str(round(MTI(ch,k)*1000)/1000),'Color','k')
                else
                    % or white number below the top of the bar
                    % (for high MTI values)
                    text(k-0.25,MTI(ch,k)-0.05, ...
                        num2str(round(MTI(ch,k)*1000)/1000),'Color',[1 1 1])
                end
            end



            subplot(3,1,3)

            % Plot of octave band signal and noise levels
            hold on

            % x-axis
            set(gca,'XTickLabel',{'125', '250', '500', '1k', '2k', '4k', '8k'})
            xlabel('Octave Band Centre Frequency (Hz)')

            % y-axis
            ymax = 10*ceil(max(Level(1,ch,:))/10);
            if ymax > 40
                ylim([0 ymax])
            end
            ylabel('SPL (dB)')

            % Sound pressure level of the received
            plot(squeeze(Level(1,ch,:)),'b','Marker','o','DisplayName','Received')
            Levels(:,1,ch) = Level(1,ch,:);

            if AuditoryMasking

                % Sound pressure level of masking from the band below
                plot(10*log10(squeeze(Iam(1,ch,:))),'Color',[0 0.5 0], ...
                    'Marker','x','DisplayName','Masking')
                Levels(:,2,ch) = 10*log10(Iam(1,ch,:));
                % Sound pressure level of the auditory reception threshold
                plot(10*log10(Irt),'k','LineStyle',':', ...
                    'DisplayName','Threshold')
                Levels(:,3,ch) = 10*log10(Irt);
                % Sum of all sources of noise
                Isum = permute(I(1,ch,:),[2,3,1]) + Irt + permute(Iam(1,ch,:),[2,3,1]);
                plot(10*log10(squeeze(Isum)), 'r', ...
                    'LineWidth',1,'LineStyle','--', ...
                    'Marker','+', 'DisplayName','Total S+N')
                Levels(:,4,ch) = 10*log10(squeeze(Isum));
            end
            legend('show','Location','EastOutside');
            hold off
        end
        
        if isstruct(in)
            mf = [0.63,0.8,1,1.25,1.6,2,2.5,3.15,4,5,6.3,8,10,12.5];
            doresultleaf(MTF_out(1:7,1:14,:),'Sparse modulation TF',{'Modulation_frequency'},...
                         'Modulation_frequency', num2cell(mf),                                'Hz',          true,...
                         'Frequency',            num2cell([125,250,500,1000,2000,4000,8000]), 'Hz',          true,...
                         'Channel',              chanID,                                      'categorical', [],...
                         'name','Modulation_TF_STI');

            if AuditoryMasking > 0
                doresultleaf(Levels,'SPL [dB]',{'Frequency'},...
                             'Frequency', num2cell([125,250,500,1000,2000,4000,8000]),                       'Hz',          true,...
                             'Level',     {'Received SPL','Masking','Threshold','Total S+N'}, 'categorical', [],...
                             'Channel',   chanID,                                                            'categorical', [],...
                             'name','Band_SPL_STI');
            else
                doresultleaf(Levels,'SPL [dB]',{'Frequency'},...
                             'Frequency', num2cell([125,250,500,1000,2000,4000,8000]), 'Hz',          true,...
                             'Level',     {'Received SPL'},                            'categorical', [],...
                             'Channel',   chanID,                                      'categorical', [],...
                             'name','Band_SPL_STI');
            end
        end

    end
else
    out = [];
end
end
%eof

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% FUNCTION TO CALCULATE MTF
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [MTF, I, Level] = STIPA_direct_MTF(audio,fs,fc,Fm)
% This is a separate function so that it can be used for the reference 
% signal (if desired) and the received signal.





% -------------------------------------------------------------------------
% FILTER INTO OCTAVE BANDS AND SQUARE
% -------------------------------------------------------------------------

% Use AARAE's linear phase octave band filters
[len, chans]= size(audio);
filterorder = 24; % very steep filter skirts improve the performance
Intensity = octbandfilter_viaFFT(audio,fs,fc,filterorder) .^2;
I = mean(Intensity); % mean square intensity
Level = 10*log10(I); % intensity expressed in decibels





% -------------------------------------------------------------------------
% SMOOTH THE ENVELOPE AND TRUNCATE
% -------------------------------------------------------------------------

% low-pass filter at 100 Hz to filter the envelope (as per standard)
% this filter also is zero phase to minimise impact on MDR
[bl, al] = butter(1,100/(0.5*fs),'low');
Intensity = filtfilt(bl, al, Intensity);

% truncate start and end to remove filter build-up period
Intensity = Intensity(round(fs/40): end-round(fs/40),:,:);
len = length(Intensity);





% -------------------------------------------------------------------------
% CALCULATE MODULATION DEPTH RATIO, MDR, and MODULATION TRANSFER FUNCTION
% -------------------------------------------------------------------------

% number of whole number cycles to use for each modulation frequency
Fm_cycles = floor(len .* Fm./fs);

% number of samples to use for each modulation frequency
Fm_len = floor(fs.*Fm_cycles ./ Fm);

% Time in seconds for each sample
t=((1:len)-1)'./fs;

% time replicated across channels
tch = repmat(t,1,chans);

mdr = zeros(1,chans,7,2);
for k = 1:7
    for j = 1:2
        % modulation depth of the received signal
        mdr(1,:,k,j) = 2 * ((sum(Intensity(1:Fm_len(k,j),:,k) ...
            .* sin(2*pi*Fm(k,j)*tch(1:Fm_len(k,j),:)))).^2 ...
            + (sum(Intensity(1:Fm_len(k,j),:,k) ...
            .* cos(2*pi*Fm(k,j)*tch(1:Fm_len(k,j),:)))).^2).^0.5 ...
            ./ sum(Intensity(1:Fm_len(k,j),:,k));
    end
    
end


m = 0.55; % Modulation depth of STIPA signal modulators

% Transform mdr to mtf using modulation depth of STIPA signal
MTF = mdr ./ m;

end