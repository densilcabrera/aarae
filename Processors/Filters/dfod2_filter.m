function [OUT, varargout] = dfod2_filter(IN, n, r,fs)
% This function call's Ivo Petras' digital fractional order differentiator
% and integrator (dfod2) code, which generates a filter. This filter is
% then applied to the input audio, which is returned to the output.
%
% Help on defod2:
% sysdfod=dfod2(n,T,r): digital fractional order differentiator
%                       and integrator    
%
% Output: =>
% Discrete system in the form of the FIR filter of the order n
% obtained by power series expansion of the backward difference.
%
% Inputs: <=
% n: order of truncation (min n=100 is recommended)
% T: sampling period in [sec]
% r: approximated fractional order (s^r), r is generally real number
%
% Author: Ivo Petras (ivo.petras@tuke.sk)
%
% Note: differentiator  -> nonrecusrive approximation 
%       integrator      -> recursive approximation
%
% Copyright (c), 2003-2011.
%

if nargin ==1 
    
    param = inputdlg({'Order of trunctation (minimum of 100 is recommended)';... 
        'Approximated fractional order (+ve for differentiation, -ve for integration)'},...
        'Dfod2 Filter',... 
        [1 60],... 
        {'100';'1'}); 
    
    param = str2num(char(param)); 
    
    if length(param) < 2, param = []; end 
    if ~isempty(param) 
        n = round(param(1));
        r = param(2);
    else
        OUT = [];
        return
    end
else
    param = [];
end
if isstruct(IN) 
    audio = IN.audio; 
    fs = IN.fs;       
elseif ~isempty(param) || nargin > 1    
    audio = IN;
end

if r == 0
    OUT = [];
    return
end

if ~isempty(audio) && ~isempty(fs) && ~isempty(n) && ~isempty(r)
    
    sysdfod = dfod2(n,1/fs,r); % make filter
    audio = filter(sysdfod.num{1},sysdfod.den{1},audio); % apply filter
    
    if isstruct(IN)
        OUT = IN; 
        OUT.audio = audio; 
        OUT.funcallback.name = 'dfod2_filter.m'; 
        OUT.funcallback.inarg = {n,r}; 
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
%  * Neither the name of the University of Sydeny nor the names of its contributors
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