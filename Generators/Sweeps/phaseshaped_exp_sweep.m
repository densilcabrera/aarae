function OUT = phaseshaped_exp_sweep(dur,start_freq,end_freq,overshoot,order,fs,reverse,phase,tukeywinratio)
% This function generates an exponential sweep and its inverse filter,
% providing for maximum phase and minimum phase solutions (as well as the
% usual zero phase solution). This may be useful in situations where the
% temporal smear of the impulse response is of particular concern (e.g. low
% frequency measurement of highly damped systems).
%
% The default parameters for this sweep are chosen for low frequency
% measurements (between 20 Hz and 500 Hz), using maximum phase. The
% resulting sweep and inverse sweep yield a very high numerical signal to
% noise ratio in the time period immediately following the impulse.
% In other words, the numerical noise is shifted away from the part of the
% impulse response that is of interest in the analysis of a system's decay
% rate (e.g. for low frequency reverberation time measurement).
%
% This is a preliminary version of the concept, and the implementation
% probably can be significantly improved. In the current method, a sweep is
% generated, and filtered using a high-order bandpass filter. The inverse
% sweep is filtered separately, so that their phase responses reinforce
% (rather than cancel) each other (in the case of minimum or maximum phase
% filters). Due to the signal-to-noise ratio limit of the bandpass filter
% (which uses cepstrum processing), the end of sweep (and inverse sweep) do
% not continue their fade without further treatment. This further treatment
% extrapolates the fade-out of the sweep and inverse sweep, which provides
% a significant improvement in signal-to-noise ratio.
%
% INPUTS
%
% Duration - in seconds, is the total duration of the waveform including
% the overshoot periods
%
% Start frequency - in Hz, is the low pass-band limit frequency of the
% filter.
%
% End frequency - in Hz, is the high pass-band limit frequency of the
% filter.
%
% Overshoot - is the number of octaves (or fractional octaves) used for the
% sweep beyond the passband. Note that this limits the possible upper
% frequency of the sweep (e.g., if the overshoot is 2 octaves, then the
% highest possible passband frequency is fs/8).
%
% Filter order - is the filter 'order' used for the filter (where the
% skirts are order*6 dB/octave). Note that since the filter is applied to
% both the sweep and inverse sweep, the skirt slope of them convolved is
% double that implied by the filter order.
%
% fs is the sampling rate in Hz.
%
% The sweep may be ascending or descending (reverse).
%
% The filter may be maximum, minimum or zero phase. Probably there is no
% reason to use this function for zero phase response (because that is the
% result from conventional sweep generation).
%
% Code by Densil Cabrera, version 1




if nargin == 0
    param = inputdlg({'Duration [s]';...
        'Start frequency [Hz]';...
        'End frequency [Hz]';...
        'Low frequency overshoot (octaves)';...
        'High frequency overshoot (octaves)';...
        'Filter order in-band';...
        'Filter order out-of-band';...
        'Tukey window fraction [0:1]';...
        'Sampling Frequency [samples/s]';...
        'Ascending [0] or descending [1] sweep';...
        'Maximum [-1], zero [0] or minimum [1] phase'},...
        'Sine sweep input parameters',1,{'60';'20';'500';'2';'1';'48';'48';'0.02';'48000';'0';'-1'});
    param = str2num(char(param));
    if length(param) < 11, param = []; end
    if ~isempty(param)
        dur = param(1);
        start_freq = param(2);
        end_freq = param(3);
        overshoot = [param(4) param(5)];
        order = [param(6) param(7)];
        tukeywinratio = param(8);
        fs = param(9);
        reverse = param(10);
        phase = param(11);
    end
else
    param = [];
end


if ~isempty(param) || nargin ~=0
    if ~exist('fs','var')
        fs = 48000;
    end
    if length(overshoot)==1
        overshoot = [overshoot overshoot];
    end
    if length(order) == 1
        order = [order order];
    end
    
    overshoot(overshoot<0.1) = 0.1;
    order(order<1) = 1;
    
    SI = 1/fs;
    ampl = 0.5;
    scale_inv = 1;
    maxfreq = (fs/2) / 2^overshoot(2); % maximum possible frequency, taking overshoot into account
    if end_freq > maxfreq, end_freq = maxfreq; end
    
    % generate sweep in time domain
    w1 = 2*pi*start_freq / 2^overshoot(1); w2 = 2*pi*end_freq * 2^overshoot(2);
    if w2 > pi*fs, w2=pi*fs; end % this should not be necessary
    K = (dur*w1)/(log(w2/w1));
    L = log(w2/w1)/dur;
    t = (0:round(dur/SI)-1)'*SI;
    phi = K*(exp(t*L) - 1);
    freq = K*L*exp(t*L);
    amp_env = 10.^((log10(0.5))*log2(freq/freq(1)));
    S = ampl*sin(phi);
    Sinv = flipud(S).*amp_env;
    
    S = bandpass(S, start_freq, end_freq, order, fs, phase);
    Sinv = bandpass(Sinv, start_freq, end_freq, order, fs, phase);
    
    % clean up the tails of the sweep & inverse sweep in the time domain -
    % they should keep decaying, but the banpass filter's noise floor limits
    % this. Let's remove the noise floor
    
    % First do the sweep
    % Find envelope, in dB
    Senv = 20*log10(abs(hilbert(S)));
    hicutind = find(abs(2*pi*end_freq - freq) == min(abs(2*pi*end_freq - freq)),1,'first');
    Sendupper = Senv(hicutind:end);
    % assume 30 dB decay exists
    Send30 = find(Sendupper <=Sendupper(1)-30,1,'first');
    if ~isempty(Send30)
        q = polyfit((1:Send30)', Sendupper(1:Send30),1); % linear regression
        extrapolatedenv = 10.^(((1:length(Sendupper))' * q(1) + q(2))./20);
        crossover = [(1:Send30)'./Send30; ones(length(Sendupper)-Send30,1)];
        S(hicutind:end) = S(hicutind:end).* (1-crossover) +...
            S(hicutind:end).* crossover.* ...
            extrapolatedenv./10.^(Senv(hicutind:end)./20);
    end
    
    % now the inverse sweep
    Senv = 20*log10(abs(hilbert(Sinv)));
    locutind = length(Senv) - find(abs(2*pi*start_freq - freq) == min(abs(2*pi*start_freq - freq)),1,'first');
    Sinvlower = Senv(locutind:end);
    Sinv30 = find(Sinvlower <=Sinvlower(1)-30,1,'first');
    if ~isempty(Sinv30)
        q = polyfit((1:Sinv30)', Sinvlower(1:Sinv30),1); % linear regression
        extrapolatedenv = 10.^(((1:length(Sinvlower))' * q(1) + q(2))./20);
        crossover = [(1:Sinv30)'./Sinv30; ones(length(Sinvlower)-Sinv30,1)];
        Sinv(locutind:end) = Sinv(locutind:end).* (1-crossover) +...
            Sinv(locutind:end).* crossover.* ...
            extrapolatedenv./10.^(Senv(locutind:end)./20);
    end
    
    
    % Final clean up of the ends in the time domain (make sure they start
    % and end at zero)
    win = tukeywin(length(S),tukeywinratio);
    S = S.*win;
    Sinv = Sinv.*win;
    
    sig_len = length(S);
    
    
    if scale_inv == 1
       fftS = fft(S);
       Sinvfft = fft(Sinv);
       mid_freq = (start_freq + end_freq)/2;
       index = round(mid_freq/(fs/sig_len));
       const1 = abs(conj(fftS(index))/(abs(fftS(index))^2));
       const2 = abs(Sinvfft(index));
       IRscalingfactor = const1/const2;
       % Sinv = Sinv * IRscalingfactor; % scaling factor is applied in
       % convolveaudiowithaudio2
    else
        IRscalingfactor = 1;
    end
    
    if reverse
        S = fliplr(S);
        Sinv = fliplr(Sinv);
    end
    
    OUT.audio = S;
    OUT.audio2 = Sinv;
    OUT.fs = fs;
    OUT.tag = ['Phase-shaped exp sweep' num2str(dur)];
    OUT.properties.dur = dur;
    OUT.properties.sig_len = sig_len;
    OUT.properties.freq = [start_freq, end_freq];
    OUT.properties.reverse = reverse;
    OUT.properties.overshoot = overshoot;
    OUT.properties.filter_order = order;
    OUT.properties.overshoot = phase;
    OUT.properties.IRscalingfactor = IRscalingfactor; % used by convolveaudiowithaudio2.m
    OUT.funcallback.name = 'phaseshaped_exp_sweep.m';
    OUT.funcallback.inarg = {dur,start_freq,end_freq,overshoot,order,fs,reverse,phase,tukeywinratio};
else
    OUT = [];
end

% This function is adapted from an exponential sweep function written by
% Nicolas Epain.
%**************************************************************************
% Copyright (c) 2015, Densil Cabrera
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