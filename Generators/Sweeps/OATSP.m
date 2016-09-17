% Generates a linear sweep (Optimized Aoshima Time Streteched Pulse)
% and its inverse filter for impulse response measurement,
% based on the concepts in:
% 
% Yôiti Suzuki, Futoshi Asano, Hack-Yoon Kim, Toshio Sone (1995)
% "An optimum computer-generated pulse signal suitable for the measurement 
% of very long impulse responses"
% Journal of the Acoustical Society of America 97(2):1119-1123
%
% The sweep bandwidth is always 0 Hz to the Nyquist frequency.
%
% code by Densil Cabrera & Daniel Jimenez
% version 1.01 (14 March 2014)

function OUT = OATSP(dur,mratio,fs,donorm,reverse)

if nargin == 0
    param = inputdlg({'Duration [s]';...
                      'Ratio of main sweep duration to remaining duration [between 0 and 1]';...
                      'Sampling rate [Hz]';...
                      'Normalize output [0 | 1]';...
                      'Ascending [0] or descending [1] sweep'},...
                      'OATSP input parameters',1,{'1';'0.5';'48000';'1';'0'});
    param = str2num(char(param));
    if length(param) < 5, param = []; end
    if ~isempty(param)
        dur = param(1);
        mratio = param(2);
        fs = param(3);
        donorm = param(4);
        reverse = param(5);
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
    
    N = 2*ceil(dur * fs/2);
    m = round((N*mratio)/2);
    k = (0:N/2)';
    Hlow = exp(1i * 4 * m * pi * k.^2 ./ N.^2);
    H = [Hlow;conj(flipud(Hlow(2:end)))];
    Sinv = ifft(H);
    Sinv = circshift(Sinv,-round(N/2-m));
    
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
    OUT.properties.mratio = mratio;
    OUT.properties.m = m;
    OUT.properties.N = N;
    OUT.properties.freq = [0, fs/2];
    OUT.properties.reverse = reverse;
    OUT.properties.IRscalingfactor = IRscalingfactor; % used by convolveaudiowithaudio2.m
    OUT.funcallback.name = 'OATSP.m';
    OUT.funcallback.inarg = {dur,mratio,fs,donorm,reverse};
else
    OUT = [];
end

%**************************************************************************
% Copyright (c) 2014, Densil Cabrera & Daniel Jimenez
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