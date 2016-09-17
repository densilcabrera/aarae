function [OUT varargout] = logspaced_complextone(fs, duration, f0, density, octavetype, phase, slope)
% This function generates a logarithmically spaced complex tone

if nargin == 0 

    param = inputdlg({'Audio sampling rate (Hz)';...
        'Tone duration (s)';...
        'Frequency of lowest tone component (Hz)';...
        'Number of components per octave';...
        'Just octaves (0) or IEC octaves (1)';...
        'Phase in degrees, or -1 for random phase';...
        'Spectral magnitude slope (dB/octave)'},...
                      'Settings',... 
                      [1 60],... 
                      {'48000';'1';'10';'3';'1';'-1';'0'}); % Default values

    param = str2num(char(param)); 

    if length(param) < 7, param = []; end 
    if ~isempty(param) 
        fs = param(1);
        duration = param(2);
        f0 = param(3);
        density = param(4);
        octavetype = param(5);
        phase = param(6);
        slope = param(7);
    end
else
    param = [];
end


if ~isempty(param) || nargin ~= 0
    
    fmax = fs/2;
    
    if octavetype == 1
        fstep = 3/density;
        f0index = 10*log10(f0);
        fmaxindex = 10*log10(fmax);
        findex = f0index:fstep:fmaxindex;
        f = 10.^(findex./10);
    else
        fstep = 1/density;
        f0index = log2(f0);
        fmaxindex = log2(fmax);
        findex = f0index:fstep:fmaxindex;
        f = 2.^findex;
    end
    
    audio = zeros(round(fs*duration),1);
    audio2 = audio;
    t = ((0:(length(audio))-1)./fs)';
    
    for k = 1:length(f)
        magnitude = 10.^((k/density * slope)/20);
        
        if phase >= 0
            audio = audio + magnitude * cos(2*pi*f(k).*t+pi*phase/180);
            audio2 = audio2 + 1/magnitude * cos(2*pi*f(k).*t+pi*phase/180);
        else
            phi = 2*rand*pi;
            audio = audio + magnitude * cos(2*pi*f(k).*t+phi);
            audio2 = audio2 + 1/magnitude * cos(2*pi*f(k).*t+phi);
        end
    end
    audio2 = flipud(audio2./max(abs(audio2)));
    audio = audio ./ max(abs(audio));
    
    if nargin == 0
        OUT.audio = audio; 
        OUT.audio2 = audio2;     
        OUT.fs = fs;       
        OUT.tag = 'logspaced tone';
        OUT.funcallback.name = 'logspaced_complextone.m';
        OUT.funcallback.inarg = {fs, duration, f0, density, octavetype, phase, slope};
    end
    
   
    if nargin ~= 0
        OUT = audio;
        varargout{1} = fs;
        %varargout{2} = ?;
    end
else
    
    OUT = [];
end

end % End of function

%**************************************************************************
% Copyright (c) 2013, Densil Cabrera
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