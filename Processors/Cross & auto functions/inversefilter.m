function OUT = inversefilter(in,fs,Threshold,hiF,loF,rolloffslope,zeropad,octsmooth,linsmooth,phasemode,rotateinv,windowoutput,writetoaudio2)
% This function derives an inverse filter of a waveform by
% inverting its spectrum. The function is interactive - allowing inspection
% and successive refinement of the inverse filter prior to finalising the
% output. Inverse filtering is not straightforward, and so this process of
% inspection and refinement is almost always necessary.
%
% Before running this function, it is advisable to carefully prepare (edit)
% the impulse response that you will use as the input.
%
% Currently this function only operates on a single channel (single
% vector), and if your input has more than 1 dimension than a dialog box
% allows you to select the vector that you wish to process.
%
% This function currently only runs using a dialog box (it can be called by
% other functions, but the inverse filter design is always interactive).
%
% The function scales the inverse filter to have the same rms amplitude as
% the input waveform.
%
% INPUT PARAMETERS:
%
% Threshold
% This is the dynamic range in dB over which spectrum inversion is done.
% Spectrum components outside this range are zeroed. This avoids dividing
% by very small numbers, dividing by zero, or dividing by noise. Set it
% extremely low (e.g. -9999) if you with to make it ineffective.
%
% hiF, loF
% A high-pass and/or low-pass filter is applied to the inverted spectrum.
% Spectrum inversion is unlikely to be useful near the Nyquist frequency or
% near 0 Hz. You should consider the characteristics of the input spectrum,
% the sampling rate, and the waveform duration in selecting these
% frequencies.
%
% rolloffslope
% The slope (in dB/octave) of the hi-pass and lo-pass filters (described
% above)
%
% zeropad
% The input can be zero-padded to increase spectrum resolution, or for
% other reasons. Zeros are concatenated to the end of the input waveform.
%
% octsmooth
% The magnitude of the inverse filter can be smoothed using a fractional
% octave band value. Use 1 for octave-band smoothing, 3 for 1/3-octave-band
% smoothing, etc.
%
% linsmooth
% The magnitude of the inverse filter can be smoothed on a linear frequency
% scale. This is done by using a Hann window as b coefficients in filtfilt
% (forward and backwards filtering, yielding zero phase). The bandwidth of
% the Hann window is specified as the linsmooth input argument. Note that
% the actual bandwidth is smaller than the Hann window length.
%
% phasemode
% This controls the phase of the inverse filter:
% 0 - zero or linear phase (ignoring the phase of the input)
% 1 - minimum phase (ignoring the phase of the input)
% -1 -  maximum phase (ignoring the phase of the input)
% 2 - complex inversion (compensating for the phase of the input)
%
% rotateinv
% The inverse filter can be rotated by half its length (using ifftshift)
%
% windowoutput
% The inverse filter can be windowed:
% 0 - no window function (i.e. rectangular window)
% 1 - Hann window function (fade-in and fade-out)
% 2 - half-Hann window function (fade-out)
% 3 - half-Hann window function (fade-in)
% 4 - Blackman window function (fade-in and fade-out)
% 5 - half-Blackman window function (fade-out)
% 6 - half-Blackman window function (fade-in)
%
% writeaudio2
% By default the inverse filter is written to the audio field of the
% output, and the input audio is written to the audio2 field. This allows
% quick testing of the inverse filter using AARAE's '*' (convolve audio with
% audio 2 button. 
% Alternatively, the inverse filter can be written to audio2, with the
% input writtent to audio. Or the inverse filter can be output without
% writing to audio2.
% 0 - write inverse filter to audio only
% 1 - write inverse filter to audio, and the input waveform to audio2
% 2 - write the input waveform to audio, and the inverse filter to audio2
%
% Code by Densil Cabrera
% version 1 (12 Sept 2016)




choose_from_higher_dimensions(in,1,1); % for now this function only operates on a vector
[bands,chans,dim4,dim5,dim6] = deal(1); %[~,bands,chans,dim4,dim5,dim6] = size(in.audio);


if nargin < 2
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
    

   
    %dialog box for settings
    prompt = {'Level threshold',...
        'High cutoff frequency (Hz)', ...
        'Low cutoff frequency (Hz)',...
        'Slope of response beyond cutoff (dB/octave)',...
        'zero pad (samples)',...
        'Fractional octave band smoothing (e.g. ''3'' for 1/3-octave) or ''0'' for no smoothing', ...
        'Linear smoothing bandwidth (Hz) - use ''0'' for no linear smoothing, and it is probably better to choose between fractional octave band and linear smoothing rather than using both',...
        'Phase response: zero phase [0], minimum phase [1], maximum phase [-1], complex inversion [2]',...
        'Rotate inverse filter via ifftshift [0 | 1]',...
        'Window the output: None [0], Hann [1], half-Hann fade-out [2], half-Hann fade-in [3], Blackman [4], half-Blackman fade-out [5], half-Blackman fade-in [6]',...
        'Write inverse filter to audio field with original in audio2 field [1] or vice versa [2], or return the inverse filter only [0]'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'-60', num2str(0.5*Nyquist), '80', '6', '0', '0', '0', '2', '0','0', '1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        Threshold = str2num(answer{1,1});
        hiF = str2num(answer{2,1});
        loF = str2num(answer{3,1});
        rolloffslope = str2num(answer{4,1});
        zeropad = str2num(answer{5,1});
        octsmooth = str2num(answer{6,1});
        linsmooth = str2num(answer{7,1});
        phasemode = str2num(answer{8,1});
        rotateinv = str2num(answer{9,1});
        windowoutput = str2num(answer{10,1});
        writetoaudio2 = str2num(answer{11,1});
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


if ~isempty(fs) && ~isempty(Threshold) && ~isempty(hiF) && ~isempty(loF) && ~isempty(zeropad) && ~isempty(rotateinv)
    while CheckResult == 1
        if zeropad > 0
            audio = [audio; zeros(zeropad,chans,bands,dim4,dim5,dim6)];
        end
        
        % make length even
        len = size(audio,1);
        if rem(len,2) == 1
            audio = [audio; zeros(1,chans,bands,dim4,dim5,dim6)];
            len = len+1;
        end
        invfilter = fft(audio, len);
        
        % fractional octave band smoothing of spectrum
        if octsmooth > 0
            magspectrum = octavesmoothing(abs(invfilter), octsmooth,fs);
            invfilter = magspectrum.*exp(1i*angle(invfilter));
        end
        
        % linear smoothing of spectrum
        if linsmooth > 0
            b = hann(round(linsmooth.*size(invfilter,1)/fs));
            invfilter = filtfilt(b,1,abs(invfilter)).*exp(1i.*angle(invfilter));
        end
        
        % Find threshold components
        magThreshold = 10.^(Threshold / 20);
        magnitude_max = repmat(max(abs(invfilter)), [len, 1,1,1,1,1]);
        below_threshold = magnitude_max < magnitude_max * magThreshold;
        invfilter = 1 ./ invfilter; % invert spectrum
        invfilter(below_threshold) = 0; % apply threshold
        
        f = fs.*(0:len-1)./len;
        order = rolloffslope*6;
        % lowpass filter
        if (hiF < Nyquist) && (hiF > loF)  && (hiF > 0)
            indhi = ceil(hiF / (fs / len)) + 1;
            %invfilter(hicomponent:len - hicomponent+2,:,:,:,:,:) = 0;
            mag = ones(size(invfilter));
            mag(indhi+1:len/2+1) = ...
            (f(indhi+1:len/2+1) ./ hiF).^-order;
            mag(len/2+2:end) = flipud(mag(2:len/2));
            mag(1) = 0;
            mag(len/2) = 0;
            invfilter = invfilter .* mag;
        end
        
        % hipass filter
        if (loF > 0) && (loF < hiF) && (loF < Nyquist)
            indlo = floor(loF / (fs / len)) + 1;
            %invfilter(1:locomponent,:,:,:,:,:) = 0;
            %invfilter(len-locomponent+2:len,:,:,:,:,:) = 0;
            mag = ones(size(invfilter));
            mag(2:indlo-1) = ...
            (f(2:indlo-1)./ loF ).^order;
            mag(len/2+2:end) = flipud(mag(2:len/2));
            mag(1) = 0;
            mag(len/2) = 0;
            invfilter = invfilter .* mag;
        end
        
        % get rid of NaN and inf if there are any
        invfilter(isnan(invfilter)) = 0;
        invfilter(isinf(invfilter)) = 0;
        
        
        % phasemode
        switch phasemode
            case 0 % zero phase
                invfilter = ifftshift(abs(invfilter));
            case 1 % minimum phase
                invfilter = minphasefreqdomain(abs(invfilter),120,1);
            case -1 % maximum phase
                invfilter = conj(minphasefreqdomain(abs(invfilter),120,1));
        end
        
        % Return to time domain, make sure it is real-valued
        out = real(ifft(invfilter));
        
        % Match rms of output to input
        out = out .* rms(audio)./rms(out);
        
        
        if rotateinv ==0 && (phasemode == 0 || phasemode == 2)
            out = ifftshift(out);
        end
        
        if rotateinv ==1 && (phasemode == 1 || phasemode == -1)
            out = ifftshift(out);
        end
        
        switch windowoutput
            case 1 % Hann window
                out = out .* hann(length(out));
            case 2
                w = hann(2*length(out));
                w = w(end/2+1:end);
                out = out.*w;
            case 3
                w = hann(2*length(out));
                w = w(1:len);
                out = out.*w;
            case 4
                out = out .* blackman(length(out));
            case 5
                w = blackman(2*length(out));
                w = w(end/2+1:end);
                out = out.*w;
            case 6
                w = blackman(2*length(out));
                w = w(1:len);
                out = out.*w;
        end
        
        % convolve input with inverse filter to test
        convlen = length(audio)+length(out)-1;
        testresult = ifft(fft(out,convlen).*fft(audio,convlen));
        outcolor = [0 0.7 0];
        incolor = [1 0 0];
        testcolor = [0.5 0 0.5];
        
        % Figure to check the result
        fig1 = figure('Name','Invert Spectrum Result');
        subplot1 = subplot(2,2,1);
        plot((0:length(out)-1)./fs,out./max(abs(out)),...
            'Color',outcolor);
        hold on
        plot((0:length(audio)-1)./fs,audio./max(abs(audio)),...
            'Color',incolor);
        hold on
        
        plot((0:convlen-1)./fs,testresult./max(abs(testresult)),'r',...
            'Color',testcolor);
        xlabel('Time (s)')
        ylabel('Amplitude')
       
        subplot2 = subplot(2,2,2);
        plot((0:length(out)-1)./fs,mag2db(abs(out)./max(abs(out))),...
            'Color',outcolor,'DisplayName','Inv filt');
        hold on
        plot((0:length(audio)-1)./fs,mag2db(abs(audio)./max(abs(audio))),...
            'Color',incolor,'DisplayName','Input');
        hold on
        
        plot((0:convlen-1)./fs,mag2db(abs(testresult)./max(abs(testresult))),'r',...
            'Color',testcolor,'DisplayName','Convolved');
        xlabel('Time (s)')
        ylabel('Level (dB)')
        ylim([-100 0])
        legend('show')

        subplot3 = subplot(2,2,3:4);
        f = fs*((1:round(length(invfilter)/2))-1)./length(invfilter);
         invfilterspectrum = fft(invfilter);
        invfilterspectrum = abs(invfilterspectrum(1:length(f))) ./ max(abs(invfilterspectrum(1:length(f))));
        semilogx(f,mag2db(invfilterspectrum),...
            'Color',outcolor);
        hold on
        
        audiospectrum = fft(audio);
        audiospectrum = abs(audiospectrum(1:length(f))) ./ max(abs(audiospectrum(1:length(f))));
        semilogx(f,mag2db(audiospectrum),...
            'Color',incolor);
        hold on
        
        testresultspectrum = fft(testresult);
        testresultspectrum = abs(testresultspectrum(1:length(f))) ./ max(abs(testresultspectrum(1:length(f))));
        semilogx(f,mag2db(testresultspectrum),'r',...
            'Color',testcolor);
        xlabel('Frequency (Hz)')
        ylabel('Level (dB)')
        ylim([-100 0])
        xlim([20 Nyquist])
        q = questdlg('OK, or adjust settings?','Make Inverse Filter','Adjust settings','Keep this inverse filter','Cancel','Keep this inverse filter');
        
        if strcmp(q,'Cancel')
            delete(fig1);
            OUT = [];
            return
        end
        
        if strcmp(q,'Adjust settings')
            
            %dialog box for settings
            prompt = {'Level threshold',...
                'High cutoff frequency (Hz)', ...
                'Low cutoff frequency (Hz)',...
                'Slope of response beyond cutoff (dB/octave)',...
                'zero pad (samples)',...
                'Fractional octave band smoothing (e.g. ''3'' for 1/3-octave) or ''0'' for no smoothing', ...
                'Linear smoothing bandwidth (Hz) - use ''0'' for no linear smoothing, and it is probably better to choose between fractional octave band and linear smoothing rather than using both',...
                'Phase response: zero phase [0], minimum phase [1], maximum phase [-1], complex inversion [2]',...
                'Rotate inverse filter via ifftshift [0 | 1]',...
                'Window the output: None [0], Hann [1], half-Hann fade-out [2], half-Hann fade-in [3], Blackman [4], half-Blackman fade-out [5], half-Blackman fade-in [6]',...
                'Write inverse filter to audio field with original in audio2 field [1] or vice versa [2], or return the inverse filter only [0]'};
            dlg_title = 'Settings';
            num_lines = 1;
            def = {num2str(Threshold), num2str(hiF), num2str(loF),...
                num2str(rolloffslope), num2str(zeropad), num2str(octsmooth),...
                num2str(linsmooth), num2str(phasemode), num2str(rotateinv),...
                num2str(windowoutput),num2str(writetoaudio2)};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            
            if ~isempty(answer)
                Threshold = str2num(answer{1,1});
                hiF = str2num(answer{2,1});
                loF = str2num(answer{3,1});
                rolloffslope = str2num(answer{4,1});
                zeropad = str2num(answer{5,1});
                octsmooth = str2num(answer{6,1});
                linsmooth = str2num(answer{7,1});
                phasemode = str2num(answer{8,1});
                rotateinv = str2num(answer{9,1});
                windowoutput = str2num(answer{10,1});
                writetoaudio2 = str2num(answer{11,1});
                CheckResult = 1;
                delete(fig1);
                audio = in.audio;
            else
                OUT = [];
                return
            end
        else
            delete(fig1);
            CheckResult = 0;
        end
        
        
    end % while
    
    
    
    if isstruct(in)
        OUT = in;
        if writetoaudio2 == 2
            OUT.audio2 = out;
        elseif writetoaudio2 == 1
            OUT.audio = out;
            OUT.audio2 = audio;
        else
            OUT.audio = out;
        end
        OUT.funcallback.name = 'inversefilter.m';
        OUT.funcallback.inarg = {fs,Threshold,hiF,loF,rolloffslope,zeropad,octsmooth,linsmooth,phasemode,rotateinv,windowoutput,writetoaudio2};
    else
        OUT = out;
    end
end