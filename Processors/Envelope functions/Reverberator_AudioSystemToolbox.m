function [OUT, varargout] = Reverberator_AudioSystemToolbox(IN,PreDelay,HighCutFrequency,Diffusion,DecayFactor,HighFrequencyDamping,WetDryMix,fs)
% This function uses Matlab's Audio System Toolbox to apply reverberation
% to the input audio.
%
% Refer to Matlab's help for 'reverberator' for information on this.

if nargin ==1 
    
    param = inputdlg({'Pre-delay (s)';... 
        'High cut-off frequency (Hz)';...
        'Diffusion (0-1)';...
        'Decay factor (0-1)';...
        'High frequency damping (0-1)';...
        'Wet-dry mix (0-1)'},...
        'Reverberator (Audio System Toolbox) settings',... 
        [1 60],... 
        {'0';'2000'; '0.5';'0.5';'0.0005';'0.3'}); 
    
    param = str2num(char(param)); 
    
    if length(param) < 6, param = []; end 
    
    if ~isempty(param) 
        PreDelay = param(1);
        HighCutFrequency = param(2);
        Diffusion = param(3);
        DecayFactor = param(4);
        HighFrequencyDamping = param(5);
        WetDryMix = param(6);
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
    dRR = reverberator('PreDelay',PreDelay,...
        'HighCutFrequency',HighCutFrequency,...
        'Diffusion',Diffusion,...
        'DecayFactor',DecayFactor,...
        'HighFrequencyDamping',HighFrequencyDamping,...
        'WetDryMix',WetDryMix,...
        'SampleRate',fs);
    [~,~,bands,dim4,dim5,dim6] = size(audio);
    for b = 1:bands
        for d4 = 1:dim4
            for d5 = 1:dim5
                for d6 = 1:dim6
                    audio(:,:,b,d4,d5,d6) = dRR(audio(:,:,b,d4,d5,d6));
                end
            end
        end
    end  
    if isstruct(IN)
        OUT = IN; 
        OUT.audio = audio; 
        OUT.funcallback.name = 'Reverberator_AudioSystemToolbox.m'; 
        OUT.funcallback.inarg = {PreDelay,HighCutFrequency,Diffusion,DecayFactor,HighFrequencyDamping,WetDryMix,fs}; 
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