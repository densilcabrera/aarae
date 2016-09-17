% This function generates a linear sweep (Optimized Aoshima Time Streteched Pulse)
% and its inverse filter for impulse response measurement,
% based on the concepts in:
% 
% Yôiti Suzuki, Futoshi Asano, Hack-Yoon Kim, Toshio Sone (1995)
% "An optimum computer-generated pulse signal suitable for the measurement 
% of very long impulse responses"
% Journal of the Acoustical Society of America 97(2):1119-1123
%
% Rather than a sweep from 0 Hz to the Nyquist frequency, this function
% uses a lower sampling rate to generate the sweep, which is then upsampled
% (along with its inverse filter). This allows the OATSP to be used for
% measurements focussing on the low frequency range.
%
% INPUTS
%
% * dur: duration in seconds. Short durations may be problematic if the
% upsampling factor is large.
%
% * upsamplingfactor: upsampling factor is a positive number, which
% is the ratio between the sampling rate of the audio and the sampling rate
% at which the OATSP signal was generated. For example, if the audio
% sampling rate is 48000 Hz, and the upsampling factor is 100, then the
% sweep will be generated using an initial sampling rate of 480 Hz (and
% hence will span 0 Hz - 240 Hz).
%
% * resample_n and resample_beta: these are filter design parameters used
% by Matlab's resample function (see Matlab's help for further
% information).
% 
% * mratio is the OATSP design parameter (refer to the paper by Suzuki et
% al).
%
% * fs: sampling rate of the output audio.
%
% * donorm: whether or not to normalize the wave.
%
% * reverse: whether or not to time-reverse the sweep (and its inverse
% filter) - which are in fact the time reverse of each other.
%
% code by Densil Cabrera & Daniel Jimenez
% version 1.0 (3 June 2015)

function OUT = OATSP_upsampled(dur,upsamplingfactor,resample_n,resample_beta,mratio,fs,donorm,reverse)

if nargin == 0
    param = inputdlg({'Duration [s]';...
                      'Upsampling factor (positive integer)';...
                      'Resample n (proportional to length of lowpass filter used for resampling)';...
                      'Resample beta (design parameter for the Kaiser window used for lowpass filter design)';...
                      'Ratio of main sweep duration to remaining duration [between 0 and 1]';...
                      'Sampling rate [Hz]';...
                      'Normalize output [0 | 1]';...
                      'Ascending [0] or descending [1] sweep'},...
                      'OATSP input parameters',1,{'60';'100';'10';'5';'0.5';'48000';'1';'0'});
    param = str2num(char(param));
    if length(param) < 8, param = []; end
    if ~isempty(param)
        dur = param(1);
        upsamplingfactor = param(2);
        resample_n = param(3);
        resample_beta = param(4);
        mratio = param(5);
        fs = param(6);
        donorm = param(7);
        reverse = param(8);
    end   
else
    param = [];
end



if ~isempty(param) || nargin ~=0
    if ~exist('fs','var')
       fs = 48000;
    end
    if ~exist('mratio','var')
       mratio = 0.5;
    end
    if ~exist('donorm','var')
       donorm = 1;
    end
    if ~exist('reverse','var')
       reverse = 0;
    end
    if ~exist('resample_beta','var')
       resample_beta = 5;
    end
    if ~exist('resample_n','var')
       resample_n = 10;
    end
    if ~exist('upsamplingfactor','var')
       upsamplingfactor = 100;
    end
    if rem(fs,upsamplingfactor) ~= 0
        upsamplingfactor = round(fs/upsamplingfactor);
    end
    
    fs1 = fs/upsamplingfactor;
    
    
    
    N = 2*ceil(dur * fs1/2);
    m = round((N*mratio)/2);
    k = (0:N/2)';
    Hlow = exp(1i * 4 * m * pi * k.^2 ./ N.^2);
    H = [Hlow;conj(flipud(Hlow(2:end)))];
    Sinv = ifft(H);
    Sinv = circshift(Sinv,-round(N/2-m));    
    Sinv = resample(Sinv,fs,fs1,resample_n,resample_beta);    
    S = flipud(Sinv);

    
    if donorm == 1
        IRscalingfactor = max(abs(S));
        S = S ./ IRscalingfactor;
    else
        IRscalingfactor = 1;
    end
    
    if reverse == 1
        S = flipud(S);
        Sinv = flipud(Sinv);
    end

    OUT.audio = S;
    OUT.audio2 = Sinv;
    OUT.fs = fs;
    OUT.tag = ['OATSP_',num2str(mratio),'_', num2str(dur)];
    OUT.properties.dur = dur;
    OUT.properties.upsamplingfactor = upsamplingfactor;
    OUT.properties.mratio = mratio;
    OUT.properties.m = m;
    OUT.properties.N = N;
    OUT.properties.resample_n = resample_n;
    OUT.properties.resample_beta = resample_beta;
    OUT.properties.freq = [0, fs1/2];
    OUT.properties.reverse = reverse;
    OUT.properties.IRscalingfactor = IRscalingfactor; % used by convolveaudiowithaudio2.m
    OUT.funcallback.name = 'OATSP_upsampled.m';
    OUT.funcallback.inarg = {dur,upsamplingfactor,resample_n,resample_beta,mratio,fs,donorm,reverse};
else
    OUT = [];
end

%**************************************************************************
% Copyright (c) 2014-15, Densil Cabrera & Daniel Jimenez
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