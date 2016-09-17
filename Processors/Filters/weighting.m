function [OUT, varargout] = weighting(IN, fs, weight)
% This function performs weighting in the frequency domain by
% acting on the magnitude of the input audio's spectrum.
%
% Version 0 (beta) - 11 March 2014

if nargin < 3 
    weight = 'a';
    param = inputdlg({'Weighting [a,b,c,d]'},...
        'Weighting',... % This is the dialog window title.
        [1 30],... 
        {'a'}); % preset answer for dialog.

    if length(param) < 1, param = []; end 
    if ~isempty(param) 
        weight = char(param(1));
    else
        OUT = [];
        return
    end
else
    param = [];
end
if isstruct(IN) 
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
else
    audio = IN;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

if ~isempty(audio) && ~isempty(fs) && ~isempty(weight)
    
    [len,chans,bands,dim4,dim5,dim6] = size(audio);
    fftlen = 2.^nextpow2(len);
    f = fs*(((1:fftlen/2)-1)/fftlen)';
    
    switch weight
        case 'a'
            tf = 10.^(2/20) .* 12200.^2 .* f.^4 ./ ...
            ((f.^2 + 20.6.^2) .* ...
            ((f.^2 + 107.7.^2) .* (f.^2 +737.9.^2)).^0.5 .* ...
            (f.^2 + 12200.^2));
        case 'b'
            tf = 10.^(0.17/20) .* 12200.^2 .* f.^3 ./ ...
            ((f.^2 + 20.6.^2) .* ...
            (f.^2 + 158.5.^2).^0.5 .* ...
            (f.^2 + 12200.^2));
        case 'c'
            tf = 10.^(0.06/20) .* 12200.^2 .* f.^2 ./ ...
            ((f.^2 + 20.6.^2) .* ...
            (f.^2 + 12200.^2));
        case 'd'
            tf = f / 6.8966888496476e-5 .* ...
            (...
            (((1037918.48-f.^2).^2 + 1080768.16.*f.^2) ./ ...
            ((9837328-f.^2).^2+11723776.*f.^2)) ./...
            ((f.^2 + 79919.29).*(f.^2+1345600))).^0.5;
        otherwise
           tf = 1+f*0;
    end
    
    tf = [tf;0;flipud(tf(2:end))];
    tf = repmat(tf,[1,chans,bands,dim4,dim5,dim6]);
    audio = ifft(fft(audio,fftlen) .* tf);
    audio = audio(1:len,:,:,:,:,:);
    if isstruct(IN)
        OUT = IN; % replicate the input structure for the output
        OUT.audio = audio; 
        OUT.funcallback.name = 'weighting.m';
        OUT.funcallback.inarg = {fs,weight};
    else
        OUT = audio;
    end
    varargout{1} = fs;
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