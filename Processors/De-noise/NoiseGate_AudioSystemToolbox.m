function [OUT, varargout] = NoiseGate_AudioSystemToolbox(IN,Threshold,AttackTime,ReleaseTime,HoldTime,fs)
% This function uses Matlab's Audio System Toolbox to apply noise gating
% (dynamics processing) to the input audio.
%
% Refer to Matlab's help for 'noiseGate' for information on this.
if nargin ==1 
    
    param = inputdlg({'Operation threshold (dB)';... 
        'Attack time (s)';...
        'Release time (s)';...
        'Hold Time (s)'},...
        'Noise Gate (Audio System Toolbox) settings',... 
        [1 60],... 
        {'-10';'0.05';'0.02';'0.05'}); 
    
    param = str2num(char(param)); 
    
    if length(param) < 4, param = []; end 
    
    if ~isempty(param) 
        Threshold = param(1);
        AttackTime = param(2);
        ReleaseTime = param(3);
        HoldTime = param(4);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
else
    param = [];
end
if isstruct(IN) 
    maxdim = 6;  % change '6' to the maximum number dimensions
    IN = choose_from_higher_dimensions(IN,maxdim,1);
    
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
elseif ~isempty(param) || nargin > 1
    audio = IN;
end


if ~isempty(audio) && ~isempty(fs)
    dRG = noiseGate('Threshold',Threshold,...
        'AttackTime',AttackTime,...
        'ReleaseTime',ReleaseTime,...
        'HoldTime',HoldTime,...
        'SampleRate',fs);
  [~,~,bands,dim4,dim5,dim6] = size(audio);
    for b = 1:bands
        for d4 = 1:dim4
            for d5 = 1:dim5
                for d6 = 1:dim6
                    audio(:,:,b,d4,d5,d6) = dRG(audio(:,:,b,d4,d5,d6));
                end
            end
        end
    end    
  
    if isstruct(IN)
        OUT = IN; 
        OUT.audio = audio; 
        OUT.funcallback.name = 'NoiseGate_AudioSystemToolbox.m'; 
        OUT.funcallback.inarg = {Threshold,AttackTime,ReleaseTime,HoldTime,fs}; 
    else
        OUT = audio;
    end
    varargout{1} = fs;
else
    OUT = [];
end

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