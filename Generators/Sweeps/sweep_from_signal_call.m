function [OUT, varargout] = sweep_from_signal_call(hiF, loF, Duration, smoothval, LevelLimit,tukeyparam,reverse)
% This function is used to call the function sweep_from_signal, which is a
% processor in AARAE (currently in the /Processors/Envelope Functions
% directory). While sweep_from_signal can be run directly, this calling
% function may be more convenient depending on the steps taken to generate
% the sweep.
%
% Sweep_from_signal derives a swept sinusoid from the power spectrum of the
% input audio. The swept sinusoid has approximately the same power spectrum
% as the input audio. An inverse filter is also output in audio2, so this
% can be used like a generator function to create a test signal for impulse
% response measurement.
%
% A possible application of this is to make a sound recording of steady
% state background noise, and to use this function to generate a sweep that
% matches the power spectrum of the background noise, allowing a
% measurement to be made with a constant signal-to-noise ratio as a
% function of frequency.
%
% It is worth experimenting with the input parameters to achieve a useable
% sweep.
%
%
% Code by Densil Cabrera
% version 0 (beta, 5 May 2014)


% Use a menu & dialog box to select a wav file or audio within AARAE
IN = choose_audio; % call AARAE's choose_audio function
if isempty(IN)
    OUT = [];
    return
end

if ~exist('hiF','var')
    hiF = IN.fs / 2.^1.16;
end
if ~exist('loF','var')
    loF = 50;
end

if nargin == 0
    
    param = inputdlg({'Upper cutoff frequency (Hz)';...
        'Lower cutoff frequency (Hz)';...
        'Duration (s)';
        'Level compensation limit (dB)';...
        'Spectrum smoothing value (2 or more for linear smoothing, or 1 or a fraction for fractional octave band smoothing)';...
        'Tukey window (fade-in and out) parameter (0:1)';...
        'Ascending [0] or descending [1] sweep'},...
        'Sweep from Signal',...
        [1 30],...
        {num2str(hiF);num2str(loF);'10';'50';'10';'0.05';'0'});
    
    param = str2num(char(param));
    
    if length(param) < 7, param = []; end
    if ~isempty(param)
        hiF = param(1);
        loF = param(2);
        Duration = param(3);
        LevelLimit = param(4);
        smoothval = param(5);
        tukeyparam = param(6);
        reverse = param(7);
    end
else
    param = [];
end


    if ~isempty(param) || nargin ~= 0
        
        OUT = sweep_from_signal(IN, hiF, loF, Duration, smoothval, LevelLimit,tukeyparam,reverse);
        
        
            
            OUT.tag = 'SweepfromSignal';      
            OUT.funcallback.name = 'sweep_from_signal_call.m';
            OUT.funcallback.inarg = {hiF, loF, Duration, smoothval, LevelLimit,tukeyparam,reverse};
            
       
    else
        % AARAE requires that in case that the user doesn't input enough
        % arguments to generate audio to output an empty variable.
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