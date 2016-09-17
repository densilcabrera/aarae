function OUT = invertspectrum(in,fs,Threshold,hiF,loF,zeropad,doweiner,rotateinv,doplot)
% This function derives a kind of 'inverse filter' of a waveform by
% inverting its spectrum. While this is almost certainly not a useful way
% to derive an inverse filter, it demonstrates the limitations of this
% simple and naive approach.
%
% Function settings are set via dialog box.
%
% A threshold (in dB) is applied to avoid dividing by numbers approaching
% zero. Set it extremely low (e.g. -9999) if you with to make it
% ineffective.
%
% A high-pass and/or low-pass filter can also be applied.
%
% The input wave can be zero-padded, which changes the resolution of the
% spectrum.
%
% Optionally, an operation a little similar to Weiner deconvolution can be
% used togenerate the inverse filter. However, this is not the same as
% Weiner deconvolution because it has no knowledge of the signal and noise
% levels that would ultimately be used for deconvolution. Instead the user
% can choose a noise-to-signal ratio (in dB) which is used for all
% frequencies. This is only included to allow experimentation.
% Use 'n' if you do not want to use this method.
%
% The filter can be rotated (so the ends shift to the middle, and the
% middle shifts to the ends). This is done using ifftshift().
%
% Code by Densil Cabrera & Daniel Jimenez
% version 0.1 (20 March 2014)


% INPUT ARGUMENTS
%
% The first input argument can either be a structure, containing the fields
% .audio (the waveform) and .fs (the sampling rate, or else it can be a
% vector or matrix containing just the audio waveform (in which case
% fs is specified by the second argument).

if nargin < 9, doplot = 0; end
if nargin < 8, rotateinv = 0; end
if nargin < 7, doweiner = 'n'; end
if nargin < 6, zeropad = 0; end
if nargin < 5, loF = 0; end
if nargin < 4, hiF = 48000; end
if nargin < 3
    if isstruct(in)
        audio = in.audio;
        fs = in.fs;
    else
        audio = in;
        if nargin < 2
            fs = inputdlg({'Sampling frequency [samples/s]'},...
                'Fs',1,{'48000'});
            fs = str2num(fs);
        end
    end
    Nyquist = fs/2;
    
    [~,bands,chans,dim4,dim5,dim6] = size(audio);
    
    %dialog box for settings
    prompt = {'Level threshold',...
        'High cutoff frequency (Hz)', ...
        'Low cutoff frequency (Hz)',...
        'Slope of response beyond cutoff (dB/octave)',...
        'zero pad (samples)',...
        'Fractional octave band smoothing (e.g. ''3'' for 1/3-octave) or ''0'' for no smoothing', ...
        'Linear smoothing bandwidth (Hz) - use ''0'' for no linear smoothing, and it is probably better to choose between fractional octave band and linear smoothing rather than using both',...
        'Phase response: simple inversion [0], minimum phase [1], maximum phase [-1]',...
        'Rotate inverse filter via ifftshift [0 | 1]',...
        'Write output to audio field [1] or to audio2 field [2]'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'-60', num2str(Nyquist), '0', '24', '0', '0', '0', '0', '0', '1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        Threshold = str2num(answer{1,1});
        hiF = str2num(answer{2,1});
        loF = str2num(answer{3,1});
        rolloffslope = str2num(answer{4,1});
        zeropad = str2num(answer{5,1});
        octsmooth = str2num(answer{6,1});
        linsmooth = str2num(answer{7,1});
        phasemode = answer{8,1};
        rotateinv = str2num(answer{9,1});
        writetoaudio2 = str2num(answer{10,1});
        CheckResult = 1;
    else
        OUT = [];
        return
    end
end
if isstruct(in)
    audio = in.audio;
    fs = in.fs;
else
    audio = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
            'Fs',1,{'48000'});
        fs = str2num(fs);
    end
end
Nyquist = fs/2;


if ~isempty(fs) && ~isempty(Threshold) && ~isempty(hiF) && ~isempty(loF) && ~isempty(zeropad) && ~isempty(doweiner) && ~isempty(rotateinv) && ~isempty(doplot)
    while CheckResult == 1
        if zeropad > 0
            audio = [audio; zeros(zeropad,chans,bands,dim4,dim5,dim6)];
        end
        
        %len = 2*(ceil(length(audio)/2)); % make len even
        len = size(audio,1);
        invfilter = fft(audio, len);
        if rem(len,2) == 0
            oddlen = 1;
        else
            oddlen = 0;
        end
        
        % fractional octave band smoothing of spectrum
        
        
        % linear smoothing of spectrum
        
        
        magThreshold = 10.^(Threshold / 20);
        
        magnitude_max = repmat(max(abs(invfilter)), [len, 1,1,1,1,1]);
        below_threshold = magnitude_max < magnitude_max * magThreshold;
        
        invfilter = 1 ./ invfilter;
        
        invfilter(below_threshold) = 0;
        
        % lowpass filter
        if (hiF < Nyquist) && (hiF > loF)  && (hiF > 0)
            hicomponent = ceil(hiF / (fs / len)) + 1;
            invfilter(hicomponent:len - hicomponent+2,:,:,:,:,:) = 0;
        end
        
        % hipass filter
        if (loF > 0) && (loF < hiF) && (loF < Nyquist)
            locomponent = floor(loF / (fs / len)) + 1;
            invfilter(1:locomponent,:,:,:,:,:) = 0;
            invfilter(len-locomponent+2:len,:,:,:,:,:) = 0;
        end
        
        % get rid of NaN and inf if there are any
        invfilter(isnan(invfilter)) = 0;
        invfilter(isinf(invfilter)) = 0;
        
        
        % phasemode
        switch phasemode
            case 1 % minimum phase
                invfilter = minphasefreqdomain(abs(invfilter),120,1);
            case 2 % maximum phase
                invfilter = conj(minphasefreqdomain(abs(invfilter),120,1));
        end
        
        % Return to time domain, make sure it is real-valued
        out = real(ifft(invfilter));
        
        if rotateinv ==1
            out = ifftshift(out);
        end
        
        % convolve input with inverse filter to test
        convlen = length(audio)+length(out)-1;
        testresult = ifft(fft(out,convlen).*fft(audio,convlen));
        outcolor = [0 0.7 0];
        incolor = [0.5 0.5 0.5];
        testcolor = [0.5 0 0.5];
        
        % Figure to check the result
        fig1 = figure('Name','Invert Spectrum Result');
        subplot1 = subplot(2,2,1);
        plot((0:length(out)-1)./fs,out,...
            'Color',outcolor);
        hold on
        plot((0:length(audio)-1)./fs,audio,...
            'Color',incolor);
        hold on
        
        plot((0:convlen-1)./fs,testresult,'r',...
            'Color',testcolor);
        xlabel('Time (s)')
        ylabel('Amplitude (dB)')
       
        subplot3 = subplot(2,2,3);
        plot((0:length(out)-1)./fs,mag2db(abs(out)),...
            'Color',outcolor);
        hold on
        plot((0:length(audio)-1)./fs,mag2db(abs(audio)),...
            'Color',incolor);
        hold on
        
        plot((0:convlen-1)./fs,mag2db(abs(testresult)),'r',...
            'Color',testcolor);
        xlabel('Time (s)')
        ylabel('Level (dB)')


        q = questdlg('Try again?','Make Inverse Filter','No');
        
        if strcmp(q,'Cancel')
            OUT = [];
            return
        end
        % need a close function too
        
        if strcmp(q,'Yes')
            
            %dialog box for settings
            prompt = {'Level threshold',...
                'High cutoff frequency (Hz)', ...
                'Low cutoff frequency (Hz)',...
                'Slope of response beyond cutoff (dB/octave)',...
                'zero pad (samples)',...
                'Fractional octave band smoothing (e.g. ''3'' for 1/3-octave) or ''0'' for no smoothing', ...
                'Linear smoothing bandwidth (Hz) - use ''0'' for no linear smoothing, and it is probably better to choose between fractional octave band and linear smoothing rather than using both',...
                'Phase response: simple inversion [0], minimum phase [1], maximum phase [-1]',...
                'Rotate inverse filter via ifftshift [0 | 1]',...
                'Write output to audio field [1] or to audio2 field [2]'};
            dlg_title = 'Settings';
            num_lines = 1;
            def = {num2str(Threshold), num2str(hiF), num2str(loF),...
                num2str(rolloffslope), num2str(zeropad), num2str(octsmooth),...
                num2str(linsmooth), num2str(phasemode), num2str(rotateinv), num2str(writetoaudio2)};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            
            if ~isempty(answer)
                Threshold = str2num(answer{1,1});
                hiF = str2num(answer{2,1});
                loF = str2num(answer{3,1});
                rolloffslope = str2num(answer{4,1});
                zeropad = str2num(answer{5,1});
                octsmooth = str2num(answer{6,1});
                linsmooth = str2num(answer{7,1});
                phasemode = answer{8,1};
                rotateinv = str2num(answer{9,1});
                writetoaudio2 = str2num(answer{10,1});
                CheckResult = 1;
                delete(fig1);
            else
                OUT = [];
                return
            end
        else
            CheckResult = 0;
        end
        
        
    end % while
    
    
    
    if isstruct(in)
        OUT = in;
        if writetoaudio2 == 2
            OUT.audio2 = out;
        else
            OUT.audio = out;
        end
        OUT.funcallback.name = 'invertspectrum.m';
        OUT.funcallback.inarg = {fs,Threshold,hiF,loF,zeropad,doweiner,rotateinv,doplot};
    else
        OUT = out;
    end
end