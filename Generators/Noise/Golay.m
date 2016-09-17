function [OUT, varargout] = Golay(fs,N,silence)
% This function generates Golay complementary sequences for impulse response
% measurement. see 
% Golay MJ. Complementary series. IRE Trans. Inf. Theory. 1961, 7:82-87.
%
% IMPORTANT:
% To analyse the recorded signal, use the function "Golay_process", which
% is in the "Cross and auto functions" Processor folder in AARAE.
%
% This code is based on the tutorial by Edgar Berdahl and Julius
% Smith "Impulse response meausrement using Golay complementary sequences",
% http://cnx.org/content/m15947/latest/ (accessed 13 February 2014).
% Adapted by Densil Cabrera for AARAE
% version 1.00 (13 February 2014)


if nargin == 0 
    
    param = inputdlg({'Audio sampling rate (Hz)';... 
        'Power of 2 for each signal length';...
        'Silence duration between signal pair (s)'},...
        'Settings',... 
        [1 60],... 
        {'48000';'15';'1'}); 
    
    param = str2num(char(param)); 
    
    if length(param) < 3, param = []; end 
    if ~isempty(param) 
        fs = param(1);
        N = param(2);
        silence = param(3);
    end
else
    param = [];
end



    

    if ~isempty(param) || nargin ~= 0
        N=abs(round(N));
        
        % remove the following if you want really long sequences!
        if N > 25
            N = 25;
        end
        N0 = N;
        
        % This part of the code is adapted from Edgar Berdahl's function
        a = [1; 1];
        b = [1; -1];

        while (N>1)
            olda = a;
            oldb = b;
            a = [olda; oldb];
            b = [olda; -oldb];
            N = N - 1;
        end
        
        
        
        gaplength = round(silence*fs);
        audio = [a;zeros(gaplength,1);b];

        if nargin == 0
            OUT.audio = audio; 
            OUT.fs = fs;
            OUT.audio2 = audio;
            OUT.tag = ['Golay',num2str(N0)];
            OUT.properties.GolayN = N;
            OUT.properties.Golaysilence = silence;
            OUT.funcallback.name = 'Golay.m';
            OUT.funcallback.inarg = {fs,N,silence};
        end
        

        if nargin ~= 0
            OUT = audio;
            varargout{1} = fs;
            varargout{2} = gaplength;
        end
    else

        OUT = [];
    end
    
end % End of function

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