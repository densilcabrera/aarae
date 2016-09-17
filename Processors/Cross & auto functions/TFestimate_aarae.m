function OUT = TFestimate_aarae(in,fs,RefChan,Threshold,Duration,doplot)
% This function calculates the transfer function between pairs of audio
% channels, using Matlab's tfestimate. Alternatively, the reference audio
% can be selected from AARAE or from file.
%
% Rather than deriving the transfer function from a single window, this
% function derives it from a series of overlapping windows, and so the
% time-variance or noise in the transfer function will cause differences
% between the results in each window (which is represented in terms of a
% reduction in coherence).
% See Matlab's help on tfestimate for information on that function. The
% transfer function can be thresholded using coherence (calculated using
% Matlab's mscohere).
%
% The output is a time response (derived via inverse Fourier transform,
% after the spectrum above the Nyquist frequency has been constructed).
%
% The 'doplot' input argument is not used in the current version.
%
% Code by Densil Cabrera
% version 0 - beta (12 August 2014)

if nargin < 6, doplot = 0; end
if nargin < 5, Duration = 2; end
if nargin < 4, Threshold = -60; end
if nargin < 3
    RefChan = 1;
    if isstruct(in)
        S = in.audio;
        fs = in.fs;
        if size(S,3)>1
            if isfield(in,'bandID')
                bandID = in.bandID;
            else
                bandID = [];
            end
        end
    else
        S = in;
        if nargin < 2
            fs = inputdlg({'Sampling frequency [samples/s]'},...
                'Fs',1,{'48000'});
            fs = str2num(char(fs));
        end
        bandID = [];
    end
    [len,chans,bands,dim4,dim5,dim6] = size(S);
    
    
    maxDuration = floor(0.25*(len) / in.fs);
    defaultDuration = 2;
    if defaultDuration > maxDuration, defaultDuration = maxDuration; end
    
    %dialog box for settings
    prompt = {'Reference channel, OR 0 for additional audio input', 'Coherence threshold [0:1]', 'Output duration (s) OR 0 for 1/4 of original duration', 'Plot (0 | 1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'1','0.9', num2str(defaultDuration),'0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        RefChan = str2num(answer{1,1});
        Threshold = str2num(answer{2,1});
        Duration = str2num(answer{3,1});
        doplot = str2num(answer{4,1});
    else
        OUT = [];
        return
    end
end

if RefChan == 0
    selection = choose_audio; % call AARAE's choose_audio function
    if ~isempty(selection)
        refaudio = selection.audio; % additional audio data
        fs2 = selection.fs; % sampling rate
        
        if ~(fs2 == fs)
            % match sampling rates
            refaudio = resample(refaudio,fs,fs2);
            disp('Sampling rate mis-match - reference audio has been re-sampled.')
        end
        [len2, chans2, bands2,dim42,dim52,dim62] = size(refaudio); % new wave dimensions
        %fftlen = max([len len2]);
    else
        OUT = [];
        return
    end
    if bands2 > 1
        if bands2 ~= bands
            refaudio = sum(refaudio,3);
            disp('multiband reference audio has been mixed down')
        else
            if isfield(selection,'bandID') && ~isempty(bandID)
                bandID2 = selection.bandID;
                if bandID2 ~= bandID
                    refaudio = sum(refaudio,3);
                    disp('multiband reference audio has been mixed down')
                end
            end
        end
    end
    if (chans2 > 1) %&& chans2 ~= chans
        prompt = {'Reference channel'};
        dlg_title = 'Reference audio settings';
        num_lines = 1;
        def = {'1'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if ~isempty(answer)
            RefChan2 = str2num(answer{1,1});
            if RefChan2 > chans2
                disp('Invalid reference channel entered - Channel 1 will be used instead')
                RefChan2 = 1;
            end
            refaudio = refaudio(:,RefChan2,:,1,1,1); % Need to select in all available dimensions
        else
            OUT = [];
            return
        end
    end
else
    if RefChan > chans
        disp('Invalid reference channel entered - Channel 1 will be used instead')
        RefChan = 1;
    end
    refaudio = S(:,RefChan,:); % need to modify!
    %fftlen = len;
end

% truncate audio to match length
if length(refaudio) > len
    refaudio = refaudio(1:len);
elseif length(refaudio) < len
    S = S(1:length(refaudio),:,:);
    len = size(S,1);
end



if Duration == 0
    outlen = round(len/4);
else
    outlen = floor(Duration * fs)+1;
    if outlen > round(len/4), outlen = round(len/4); end
end




%Threshold = 0.9; % coherence threshold
% The window length should be at least twice the deired impulse
% response length
winlength = 2^nextpow2(outlen * 2)-1;

[TF,COH] = deal(zeros(ceil(winlength/2)+1,chans,bands));

% Window overlap in samples
overlap = ceil(winlength*0.9); % 90% overlap, in samples

for ch = 1:chans
    for b = 1:bands
        
        % TF includes components from DC to the Nyquist frequency
        TF(:,ch,b) = tfestimate(refaudio,S(:,ch,b),hann(winlength),overlap);
        
        % Calculate coherence function for thresholding.
        COH(:,ch,b) = mscohere(refaudio,S(:,ch,b),hann(winlength),overlap);
        
    end
end


TF(COH < Threshold) = 0; % zero all values below coherence threshold

% Include values above Nyquist frequency for the ifft.
TF = [TF;conj(flipdim(TF(2:end-1,:,:,:,:,:),1))];


% Return to time domain
out = ifft(TF);

% Truncate to the desired length
out = out(1:outlen,:,:,:,:,:);

if isstruct(in)
    OUT = in; % replicate input structure
    OUT.audio = out;
    OUT.funcallback.name = 'TFestimate_aarae.m';
    OUT.funcallback.inarg = {fs,RefChan,Threshold,Duration,doplot};
else
    OUT = out;
end

% these plots are unnecessary within AARAE, and probably need tweaking for
% them to work (to do)
if false
%if doplot
    figure('Name', 'Transfer Function')
    
    subplot(3,1,1)
    t = ((1:outlen)'-1)/in.fs;
    plot(repmat(t, [1,chans,bands]),out)
    xlabel('Time (s)')
    xlim([0 t(end)])
    ylabel('Amplitude')
    
    subplot(3,1,2)
    plot(repmat(t, [1,chans,bands]),10*log10(out.^2))
    xlabel('Time (s)')
    xlim([0 t(end)])
    ylabel('Level (dB)')
    
    subplot(3,1,3)
    lenTF = length(TF);
    TFhalf = TF(1:ceil(lenTF/2),:,:);
    lenTF_half = length(TFhalf);
    f = in.fs .* ((1:lenTF_half)'-1) ./ lenTF;
    semilogx(repmat(f(2:end),[1,chans, bands]),10*log10(abs(TFhalf(2:end,:,:)).^2))
    xlabel('Frequency (Hz)')
    xlim([20 f(end)])
    ylabel('Gain (dB)')
    
    figure
    semilogx(f(2:end),COH(3:end),'r')
    xlabel('Frequency (Hz)')
    xlim([f(2) f(end)])
    ylim([0 1])
    title('Coherence')
end


