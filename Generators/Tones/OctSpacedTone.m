function OUT = OctSpacedTone(duration,frequency,fs,tukeyratio,hannexp)
% Generates an octave spaced tone, i.e. a Shepard tone
%
% Octave spaced tones are mainly used for demonstrations of octave
% similarity vs chroma, which is exemplified by the 'tritone paradox'.


if nargin == 0
    param = inputdlg({'Duration [s]';...
                       'Frequency [Hz]';...
                       'Sampling Frequency [Hz]';...
                       'Fade-in, fade-out Tukey window ratio (0 - 1)';...
                       'Hann window exponent for spectrum envelope'},...
                       'Sine tone input parameters',1,...
                       {'10';'440';'48000';'0.1';'1'});
    param = str2num(char(param));
    if length(param) < 5, param = []; end
    if ~isempty(param)
        duration = param(1);
        frequency = param(2);
        fs = param(3);
        tukeyratio = param(4);
        hannexp = param(5);
    end
end
if ~isempty(param) || nargin ~= 0
    f = 2^((log2(frequency)-floor(log2(frequency)))+4);
    f = repmat(f,[1,11]).*2.^(0:10);
    f = f(f<fs/2);
    f = f(f<24000);
    t = linspace(0,duration,fs*duration);
    y = zeros(length(t),1);
    env = (1-0.5*(cos(2*pi*(log2(f./440)-log2(16/440))./10.5507)+1)).^hannexp;
    env = env ./ sum(env);
    for n = 1:length(f)
        y = y + cos(2*pi*f(n).*t') .* env(n);
    end
    y = y .* tukeywin(length(y),tukeyratio);
    tag = [num2str(frequency) 'Hz OctSpacedTone'];
    
    OUT.audio = y;
    OUT.fs = fs;
    OUT.tag = tag;
    OUT.funcallback.name = 'OctSpacedTone.m';
    OUT.funcallback.inarg = {duration,frequency,fs,tukeyratio,hannexp};
else
    OUT = [];
end

end % End of function

%**************************************************************************
% Copyright (c) 2016, Densil Cabrera
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