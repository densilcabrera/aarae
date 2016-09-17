function [OUT, varargout] = infrasound_octbandfilter(IN, fs, param)
% This function does octave band filtering in the infrasound range.
% It calls AARAE's octbandfilter to do the filtering, after resampling the
% audio to 192 Hz. By sending a fake sampling rate of 192 kHz to
% octbandfilter, the result is octave bands centred on frequencies between
% 0.0315 Hz to 16 Hz.
%
% It may be worthwhile to investigate the optimum length of the
% anti-aliasing filter used in downsampling.
%
% Code by Densil Cabrera
% Version 1.01 (29 December 2013)


% -------------------------------------------------------------------------
if isstruct(IN) 
    audio = IN.audio; 
    fs = IN.fs;       
    ok = 1;
    % list dialog
    if nargin < 3
        nominalfreq = [31.5,63,125,250,500,1000,2000,4000,8000,16000] ./ 1000;
        param = nominalfreq;
        [S,ok] = listdlg('Name','Octave band filter input parameters',...
                         'PromptString','Center frequencies [Hz]',...
                         'ListString',[num2str(param') repmat(' Hz',length(param),1)]);
        param = param(S);
    end
    
% ---------------------------------------------------------------------    
elseif nargin > 1
    
    audio = IN;
    param = [31.5,63,125,250,500,1000,2000,4000,8000,16000] ./ 1000;
    ok = 1;
end





if ~isempty(audio) && ~isempty(fs) && ok==1

    bands = size(audio,3);
    
    % mixdown multiband audio
    if bands > 1
        audio = sum(audio,3);
        disp('Multiband audio has been mixed into a single band')
    end
    
    
     % downsample the audio to fs = 192 Hz
    
    n = 100; % controls the length of anti-aliasing FIR filter
    audio = resample(audio,19200,fs,n);
    audio = resample(audio,1920,19200,n);
    audio = resample(audio,192,1920,n);
    
    % call octbandfilter, using a fake sampling rate and centre frequencies
        audio = octbandfilter(audio,192000,...
        param*1000);

    
    if isstruct(IN)
        OUT = IN; % replicate the input structure for output
        OUT.audio = audio; 
        OUT.fs = 192; % 192 Hz sampling rate
        OUT.bandID = param;
        OUT.funcallback.name = 'infrasound_octbandfilter.m';
        OUT.funcallback.inarg = {fs, param};
    else
        
        OUT = audio;
    end
    varargout{1} = fs;
    
else
    
    OUT = [];
end

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