function [OUT, varargout] = sweep_from_signal(IN, hiF, loF, Duration, smoothval, LevelLimit,tukeyparam,reverse,fs)
% This function derives a swept sinusoid from the power spectrum of the
% input audio. The swept sinusoid has approximately the same power spectrum
% as the input audio. An inverse filter is also output in audio2, so this
% can be used like a generator function to create a test signal for impulse
% response measurement.
%
% A possible application of this is to make a sound recording of steady
% state background noise, and to use this function to generate a sweep that
% matches the power spectrum of the background noise, allowing a
% measurement to be made with a constant signal-to-noise ratio as a
% function of frequency.
%
% It is worth experimenting with the input parameters to achieve a useable
% sweep.
%
%
% Code by Densil Cabrera
% version 0 (beta, 5 May 2014)




if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,1,1); % only 1-D data allowed
    signal = IN.audio;
    fs = IN.fs;
else
    signal = IN;
end

if ~exist('hiF','var')
    hiF = fs / 2.^1.16;
end
if ~exist('loF','var')
    loF = 50;
end


if nargin ==1
    
    param = inputdlg({'Upper cutoff frequency (Hz)';...
        'Lower cutoff frequency (Hz)';...
        'Duration (s)';
        'Level compensation limit (dB)';...
        'Spectrum smoothing value (2 or more for linear smoothing, or 1 or a fraction for fractional octave band smoothing)';...
        'Tukey window (fade-in and out) parameter (0:1)';...
        'Ascending [0] or descending [1] sweep'},...
        'Sweep from Signal',...
        [1 30],...
        {num2str(hiF);num2str(loF);'10';'50';'10';'0.05';'0'});
    
    param = str2num(char(param));
    
    if length(param) < 7, param = []; end
    if ~isempty(param)
        hiF = param(1);
        loF = param(2);
        Duration = param(3);
        LevelLimit = param(4);
        smoothval = param(5);
        tukeyparam = param(6);
        reverse = param(7);
    end
else
    param = [];
end


if ~isempty(signal) && ~isempty(fs)
    
    [len,chans,bands] = size(signal);
    
    
        if bands > 1
            signal = sum(signal,3);
            disp('Multiband audio has been mixed into a single band')
        end
    
        if chans > 1
            signal = mean(signal,2);
            disp('Multichannel audio has been mixed')
        end
    
    if ~isreal(signal)
        signal = abs(signal).*cos(angle(signal));
        disp('The inverse Hilbert transform of the complex signal is being used')
    end
    
    if hiF > fs / 2.^1.16;
        hiF = fs / 2.^1.16;
    end
    if (loF > hiF) || (loF < 10)
        loF = 50;
    end
    
    %     minfftlen = 4*fs/loF;
    %     if len < minfftlen
    %         fftlen = minfftlen;
    %     else
    %         fftlen = len;
    %     end
    
    winlen = Duration * fs * 2;
    if winlen > len
        signal = [signal;zeros(winlen-len,1)];
    end
    % use Matlab's spectrogram to get the spectrum or average of short term
    % spectra
    [spectrum, freq] = spectrogram(signal,hann(winlen),...
        round(winlen/2),winlen,fs);
    powspectrum = mean(abs(spectrum).^2,2);
    clear spectrum
    
    
    
    %SMOOTH SPECTRUM (BEFORE TRUNCATING)
    if smoothval > 1
        % Spectrum smoothing using a linear scale
        b = hann(round(smoothval)+2).^0.5; % smoothing filter (FIR)
        b = b(2:end-1); % get rid of zeros
        b = b ./sum(b); % normalize for constant gain
        powspectrum = filtfilt(b,1,powspectrum); % zero phase filter
    elseif (smoothval > 0) && (smoothval <= 1)
        % fractional octave band smoothing
        % NEED A FASTER METHOD!!
        % Probably do logarithmic resampling (reduce the number of spectrum
        % components at high frequencies) - but will need to change
        % subeseqent code to allow for non-linear sampling
        % upper cut-off freq for each band
        freqhi = freq .* 10.^(0.3*smoothval/2);
        % lower cut-off freq for each band
        freqlo = freq ./ 10.^(0.3*smoothval/2);
        
        for k = 2:length(freq)
            powspectrum(k) = mean(powspectrum(freq >= freqlo(k) & freq <= freqhi(k)));
        end
    end
    
    
    % discard spectrum outside of sweep limits (and above Nyquist)
    powspectrum = powspectrum(freq >= loF & freq <= hiF);
    freq = freq(freq >= loF & freq <= hiF);
    fstep = freq(2)-freq(1); % should be equal to inverse of window period
    len = length(powspectrum);
    
    % limit spectrum range and divide by minimum
    maxps = max(powspectrum);
    rangeps = 10.^(LevelLimit/10);
    powspectrum(powspectrum < maxps/rangeps) = maxps/rangeps;
    powspectrum = powspectrum ./ min(powspectrum);
    
    instf = zeros(ceil(sum(powspectrum)),1);
    findex = 1;
    lenk = diff([0;round(cumsum(powspectrum))]);
    for k = 1:len
        instf(findex:findex+lenk(k)-1) = freq(k) ...
            + (((1:lenk(k))'-0.5 - 0.5*lenk(k))./lenk(k)) * fstep;
        findex = findex + lenk(k);
    end
    instf = instf(1:findex-1);
    %figure; plot(instf);
    instf = resample(instf,round((20+Duration.*fs)/1000),round(length(instf)/1000),10);
    instf = instf(10:end-10);
    
    
    
    %audio = cos(cumsum(instf/(fs/2/pi)));
    audio = sin(cumsum(instf/(fs/2/pi)));
    audiospectrum = fft(audio);
    audiospectrum = conj(audiospectrum ./ abs(audiospectrum).^2);
    audiospectrumf = fs.*((1:length(audio))-1) ./ length(audio);
    audiospectrum(audiospectrumf<loF) = 0;
    audiospectrum(audiospectrumf > fs - loF) = 0;
    audiospectrum((audiospectrumf > hiF) & (audiospectrumf < fs - hiF)) = 0;
    audio2 = ifft(audiospectrum);
    
    audio = audio .* tukeywin(length(audio),tukeyparam);
    %audio2 = audio2 .* tukeywin(length(audio),tukeyparam); %probably unnecessary
    if reverse == 1
        audio = flipud(audio);
        audio2 = flipud(audio2);
    end
    
    if isstruct(IN)
        OUT = IN;
        OUT.audio = audio;
        OUT.audio2 = audio2;
        OUT.properties.dur = Duration;
        OUT.properties.LevelRange = LevelLimit;
        OUT.properties.tukey = tukeyparam;
        OUT.properties.smoothing = smoothval;
        OUT.properties.freq = [loF, hiF];
        OUT.properties.reverse = reverse;
        OUT.funcallback.name = 'sweep_from_signal.m';
        OUT.funcallback.inarg = {hiF, loF, Duration, smoothval, LevelLimit,tukeyparam,reverse};
    else
        OUT = audio;
    end
    varargout{1} = fs;
    varargout{2} = audio2;
else
    OUT = [];
end

%**************************************************************************
% Copyright (c) 2014, Densil Cabrera
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