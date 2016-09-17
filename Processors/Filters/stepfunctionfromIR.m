function OUT = stepfunctionfromIR(IN,dur,method,fs)
% This function filters the input audio with ones, matching the length of
% the input. Hence it may be used to derive the step function from an
% impulse response. 
%
% INPUTS:
% dur - output duration (and length of the ones vector) in seconds
%
% method - either filter in the time domain (method 1) or multiply in the
% frequency domain (method 2). Method 2 is usually faster (but uses more
% memory)
%
% fs -  audio sampling rate


if isstruct(IN) 
    audio = IN.audio;
    fs = IN.fs;
elseif  nargin > 1
    audio = IN;
end
[len,chans,bands,dim4,dim5,dim6] = size(audio);
if nargin ==1 
    param = inputdlg({'Output duration';... 
        'Method: Filter with ones [1], Frequency domain multiplication [2]'},...
        'Step Function from IR Settings',...
        [1 60],...
        {num2str(len./fs);'2'});
    param = str2num(char(param)); 
    if length(param) < 2, param = []; end 
    if ~isempty(param) 
        dur = param(1);
        method = param(2);
    else
        OUT = [];
        return
    end
end


if ~isempty(audio) && ~isempty(method)
    outlen = round(dur*fs);
    if outlen < len
        audio = audio(1:outlen,:,:,:,:,:);
    elseif outlen > len
        audio = [audio;zeros(outlen-len,chans,bands,dim4,dim5,dim6)];
    end
    if method == 2
        audio = ifft(fft(audio,2*outlen) .* fft(ones(size(audio)),2*outlen));
        audio = audio(1:outlen,:,:,:,:,:);
    else
        audio = filter(ones(outlen,1),1,audio);
    end
    if isstruct(IN)
        OUT = IN; 
        OUT.audio = audio; 
        OUT.funcallback.name = 'stepfunctionfromIR.m'; 
        OUT.funcallback.inarg = {dur,method,fs}; 
    else
        OUT = audio;
    end
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