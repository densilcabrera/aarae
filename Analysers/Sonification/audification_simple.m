function out = audification_simple(in,fs,hiF,loF,speedupfactor,reverse,envelopeexp,envelopesmooth,spectrumexp,spectrumsmooth)
% Provides a simple framework for audification of input audio
%
% Audification means 'playing data as if it is sound'. When we audify an
% audio waveform, this can be as trivial as playing it. However, it is also
% possible to do some simple transformations of the audio data so that
% certain features become more audible.
%
% Common transformations that are used in audification include:
% * filtering
% * speeding up or slowing down the wave (resampling it);
% * time-reversing the waveform;
% * doing envelope manipulations;
% * doing spectrum manipulations;
%
% Code by Densil Cabrera
% version 1.01 (10 February 2014)
%
% ABOUT THE SETTINGS
% The default values result in no change in the audio
%
% A high-pass and/or low-pass filter can be applied. The specified
% frequencies refer to the original file (prior to resampling if that is
% done). Filtering is done via fft.
%
% Speed-up (resampling) factor may be used to speed up or slow down the
% audio data. A value of 2 doubles the speed (with a +1 octave frequency
% shift). A value of 0.5 halves the speed (with a -1 octave frequency
% shift).
%
% Time reversal can be set (1) or clear (0, meaning no change).
%
% The envelope dynamic contrast exponent controls the time envelope of the
% audio (via the Hilbert transform):
%   1 creates no change
%   A larger number (e.g. 1) increases dynamic contrast)
%   0 creates a constant envelope (no dynamic contrast)
%   A negative value inverts the dynamic contrast.
%
% The envelope smoothing filter length is the length (in samples) of a
% filter used to smooth the envelope function. It should be a positive
% integer. A value of 1 creates no smoothing.
%
% The spectrum dynamic contrast exponent acts similarly to the envelope
% dynamic contrast exponent, except that it is applied to the spectrum.
%   1 creates no change
%   A greater number (e.g. 2) increases dynamic contrast in the spectrum)
%   0 creates a constant magnitude spectrum (no dynamic contrast)
%   A negative value inverts the magnitude spectrum.
%
% The spectrum smoothign filter works in the same way as the envelope
% smoothing filter, except that it is applied to the spectrum magnitude.
%
% Code by Densil Cabrera
% version 1.01 (10 February 2014)


if isstruct(in)
    in = choose_from_higher_dimensions(in,3,1); 
    % required field of input structure
    out = in;
    audio = in.audio;
    audio = sum(audio,3); % audio waveform, 2-d max
    fs = in.fs; % audio sampling rate
    Nyquist  = fs/2;
else
    out = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
        Nyquist = fs/2;
    end
end

if nargin < 10, spectrumsmooth = 1; end
if nargin < 9, spectrumexp = 1; end
if nargin < 8, envelopesmooth = 1; end
if nargin < 7, envelopeexp = 1; end
if nargin < 6, reverse = 0; end
if nargin < 5, speedupfactor = 1; end
if nargin < 4, loF = 0; end
if nargin < 3
    hiF = Nyquist;
    % dialog box for settings
    prompt = {'High cutoff frequency (Hz):', ...
        'Low cutoff frequency (Hz):', ...
        'Speed-up (resampling) factor:', ...
        'Time reversal (0 | 1):', ...
        'Envelope dynamic contrast exponent:', ...
        'Envelope smoothing filter length:', ...
        'Spectrum magnitude contrast exponent:', ...
        'Spectrum smoothing filter length:'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {num2str(Nyquist),'0','1','0','1','1','1','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if isempty(answer)
        out = [];
        return
    else
        hiF = str2num(answer{1,1});
        loF = str2num(answer{2,1});
        speedupfactor = str2num(answer{3,1});
        reverse = str2num(answer{4,1});
        envelopeexp = str2num(answer{5,1});
        envelopesmooth = str2num(answer{6,1});
        spectrumexp = str2num(answer{7,1});
        spectrumsmooth = str2num(answer{8,1});
    end
    envelopesmooth = abs(round(envelopesmooth-1))+1;
    spectrumsmooth = abs(round(spectrumsmooth-1))+1;
end

if ~isempty(audio) && ~isempty(fs) && ~isempty(hiF) && ~isempty(loF) && ~isempty(speedupfactor) && ~isempty(reverse) && ~isempty(envelopeexp) && ~isempty(envelopesmooth) && ~isempty(spectrumexp) && ~isempty(spectrumsmooth)
    % filtering
    if ~(hiF == Nyquist) || ~(loF == 0);
        % derive spectrum, even length
        fftlen = 2*ceil(length(audio)/2);
        spectrum = fft(in.audio,fftlen);

        % lowpass filter
        if (hiF < Nyquist) && (hiF > loF)  && (hiF > 0)
            hicomponent = ceil(hiF / (in.fs / fftlen)) + 1;
            spectrum(hicomponent:fftlen - hicomponent+2,:) = 0;
        end

        % hipass filter
        if (loF > 0) && (loF < hiF) && (loF < Nyquist)
            locomponent = floor(loF / (in.fs / fftlen)) + 1;
            spectrum(1:locomponent,:,:) = 0;
            spectrum(fftlen-locomponent+2:fftlen,:) = 0;
        end

        % return to time domain
        audio = ifft(spectrum);
    end


    % resampling
    if ~(speedupfactor == 1)
        audio = resample(audio,fs,fs*speedupfactor);
    end

    % time reversal
    if reverse
        audio = flipdim(audio,1);
    end

    % envelope contrast and smoothing
    if ~(envelopeexp == 1) || ~(envelopesmooth == 1)
        analytic = hilbert(audio);
        envelope = abs(analytic) .^ envelopeexp;
        if ~(envelopesmooth == 1)
            b = ones(1,envelopesmooth)/envelopesmooth;  % averaging filter
            envelope = fftfilt(b,envelope);% smooth the envelope
        end
        % adjust audio by envelope
        audio = envelope .* cos(angle(analytic));
    end

    % spectrum contrast
    if ~(spectrumexp == 1) || ~(spectrumsmooth == 1)
        spectrum = fft(audio);
        magnitude = abs(spectrum).^(spectrumexp);
        phase = angle(spectrum);
        %phase = angle(spectrum);
        if ~(spectrumsmooth == 1)
            b = ones(1,spectrumsmooth)/spectrumsmooth;  % averaging filter
            magnitude = fftfilt(b,magnitude);% smooth the envelope
        end
        spectrum = magnitude .* exp(1i * phase);
        audio = real(ifft(spectrum));
    end

    %normalize
    audio = audio ./max(max(abs(audio)));
    if isstruct(in)
        out.audio = audio;
        out.funcallback.name = 'audification_simple.m';
        out.funcallback.inarg = {fs,hiF,loF,speedupfactor,reverse,envelopeexp,envelopesmooth,spectrumexp,spectrumsmooth};
    else
        out = audio;
    end
else
    out = [];
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
%  * Neither the name of The University of Sydney nor the names of its contributors
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
