function OUT = SilenceSweep_Farina2009(n,fadeindur,fadeoutdur,fs,gapdur,deltaL)
% This function generates the 'silence sweep' as described by Angelo
% Farina:
%
% A. Farina (2009) "Silence Sweep: a novel method for measuring
% electro-acoustical devices," 126th AES Convention, Munich, Germany
%
% The test signal consists of a period of silence, then the silence sweep,
% then the exponential sweep. The duration of the gap between these 3
% elements is specified as a user parameter. The default parameters are
% those used in the paper, however a more sensitive analysis can be achived
% by using a longer fade-in and fade-out duration for the sweep (e.g. 1 s,
% perhaps with a higher order MLS). However, note that the fade-in and
% fade-out durations limit the effective frequency range of the sweep. The
% sensitivity of the test signal can be explored by analysing the test
% signal itself (without playing it through an audio system).
%
% Although Farina puts the silence in the middle, in this implementation
% the silence is first because this avoids any possibility that system
% reverberation could intrude into the silence.
%
% Note that in his presentation, Angelo Farina speculated that the rms
% level of the sweep should perhaps be 8-12 dB less than that of the
% silence sweep, although in his examples he states that they were the same
% rms level (creating a strange result where 2nd order harmonic distortion
% from the sweep is greater than THD+N from the silence sweep). This can be
% controlled using the last user parameter (deltaL), which attenuates the
% sweep, regardless of its sign.
%
% Like all measurements of non-linear distortion, it is best to make the
% measurement in quiet conditions, because otherwise it may be difficult to
% distinguish distortion from noise. The analysis assumes that background
% noise is steady state, and significant time-variance or impulsiveness in
% the noise may introduce artefacts in the analysed results.
%
% Use the SilenceSweepAnalysis analyser (in Non-LTI analysis) to analyse
% recordings made with this test signal. Alternatively, you can analyse it
% 'manually' by using the '*' (convolve audio with audio2) button.
%
% This method of assessing distortion may be helpful because it allows the
% distortion due to broadband excitation to be assessed (rather than
% distortion from one or two excitation frequencies). The silence sweep
% does not distinguish between harmonic and other types of non-linear 
% distortion.
%
% Code by Densil Cabrera, December 2015




if nargin == 0
    
    param = inputdlg({'MLS order (e.g., 16, 17, or 18)';...
        'Sweep fade-in duration (s)';...
        'Sweep fade-out duration (s)';...
        'Sampling rate (Hz) - must be at least 44100';...
        'Duration of the gap between silence, silence sweep and sweep (s) - cannot be less than 1 MLS';...
        'Level difference (Leq) between MLS (silence sweep) and sweep (dB)'},...
        'Silence Sweep',... 
        [1 60],... 
        {'17';'0';'0.1';'48000';'1';'0'}); 
    
    param = str2num(char(param)); 
    
    if length(param) < 6, param = []; end 
    if ~isempty(param) 
        n = round(param(1));
        fadeindur = param(2);
        fadeoutdur = param(3);
        fs = round(param(4));
        gapdur = param(5);
        deltaL = abs(param(6));
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
end


if fs < 44100, fs = 44100; end

if ~isempty(gapdur) && ~isempty(fs)
    
    % Generate MLS repeated sequence
    %n = 17;
    cycles = 41;
    mls = GenerateMLSSequence(cycles, n, 0);
    
    mlslen = 2^n-1;
    mls(mlslen*floor(cycles/2)+1:mlslen*ceil(cycles/2)) = 0;
    if gapdur < round(mlslen/fs), gapdur = round(mlslen/fs); end
    
    % generate sweep 10 octaves
    % code adapted from prior code by Nicholas Epain
    SI = 1/fs;
    dur = 10*mlslen*SI;
    ampl = 0.25;
    start_freq = 20; % Hz
    end_freq = 20480; % Hz
    rcos_ms1 = fadeindur*1000; % fade-in duration in ms
    rcos_ms2 = fadeoutdur*1000; % fade-out duration in ms
    scale_inv = 1;
    w1 = 2*pi*start_freq; w2 = 2*pi*end_freq;
    K = (dur*w1)/(log(w2/w1));
    L = log(w2/w1)/dur;
    t = (1:10*mlslen-1)*SI;
    phi = K*(exp(t*L) - 1);
    freq = K*L*exp(t*L);
    amp_env = 10.^((log10(0.5))*log2(freq/freq(1)));
    S = ampl*sin(phi);
    rcos_len1 = round(length(S)*((rcos_ms1*1e-3)/dur));
    rcos_len2 = round(length(S)*((rcos_ms2*1e-3)/dur));
    sig_len = length(S);
    rcoswin1 = hann(2*rcos_len1).';
    rcoswin2 = hann(2*rcos_len2).';
    S = [S(1:rcos_len1).*rcoswin1(1:rcos_len1),S(rcos_len1+1:sig_len-rcos_len2),S(sig_len-rcos_len2+1:sig_len).*rcoswin2((rcos_len2+1):(rcos_len2*2))]';
    Sinv = flip(S).*amp_env';

    % correction for allpass delay
    Sinvfft = fft(Sinv);
    Sinvfft = Sinvfft.*exp(1i*2*pi*(0:(sig_len-1))*(sig_len-1)/sig_len)';
    Sinv = real(ifft(Sinvfft));

    if scale_inv == 1
       fftS = fft(S);
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
    
    % convolve mls (incl silence gap) with sweep
    fftlen = length(mls) + length(S) - 1;
    gaplen = round(fs*gapdur);
    audio = [zeros(gaplen,1); ifft(fft(mls,fftlen) .* fft(S,fftlen)); zeros(gaplen,1)];
    audio = 0.5.^0.5 * ampl * audio ./ ...
        (rms(audio(length(S)+gaplen+1:end-(length(S)+gaplen))));
    audio(1:length(S)+gaplen) = 0; % silence
    audio(1+end-(length(S)+gaplen):end) = [zeros(gaplen,1);S./db2mag(deltaL)]; % sweep
    
    OUT.fs = fs;
    OUT.audio = audio;
    OUT.audio2 = Sinv;
    OUT.tag = 'SilenceSweep';
    OUT.funcallback.name = 'SilenceSweep_Farina2009.m'; 
    OUT.funcallback.inarg = {n,fadeindur,fadeoutdur,fs,gapdur,deltaL};
    OUT.properties.IRscalingfactor = IRscalingfactor;
    OUT.properties.MLSorder = n;
    OUT.properties.gapdur = gapdur;
    OUT.properties.deltaL = deltaL;
    OUT.properties.fadeindur = fadeindur;
    OUT.properties.fadeoutdur = fadeoutdur;
    OUT.properties.start_freq = start_freq;
    OUT.properties.end_freq = end_freq;
    OUT.properties.sweepdur = dur;
    OUT.properties.SilenceSweep = 1; % currently the analyser just checks that this exists. In future this could have more meaning.
else
    OUT = [];
end


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
%  * Neither the name of the The University of Sydney nor the names of its contributors
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