function [OUT, varargout] = Gamma_tone_Fast_Hummersone2011(IN,fs,hif,lof,bands,align)
% This function filters the input audio using Christopher Hummersone's
% GammatoneFast function.

if nargin ==1 

    param = inputdlg({'High Centre Frequency (Hz)';... % These are the input box titles in the
                      'Low Centre Frequency (Hz)';...
                      'Number of Centre Frequencies';...
                      'Phase alignment (0 | 1)'},...% inputdlg window.
                      'Settings',... % This is the dialog window title.
                      [1 30],... 
                      {'16000';'50';'40';'0'}); % Defaults

    param = str2num(char(param)); % Since inputs are usually numbers it's a
                                  % good idea to turn strings into numbers.

    if length(param) < 4, param = []; end % You should check that the user 
                                          % has input all the required
                                          % fields.
    if ~isempty(param) % If they have, you can then assign the dialog's
                       % inputs to your function's input parameters.
        hif = param(1);
        lof = param(2);
        bands = param(3);
        align = param(4);
    else
        OUT = [];
        return
    end
    if align == 1, align = true; else align = false; end
else
    param = [];
end
if isstruct(IN) 
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
elseif ~isempty(param) || nargin > 1
  
    audio = IN;

end


if ~isempty(audio) && ~isempty(fs)
    
    cfs = MakeErbCFs(lof,hif,bands);
    
    audio = sum(audio,3); % mixdown 3rd dimension if it exists
    [len,chans,~,dim4,dim5,dim6] = size(audio);
    bm = zeros(len,chans,bands,dim4,dim5,dim6);
    for ch = 1:chans
        for d4 = 1:dim4
            for d5 = 1:dim5
                for d6 = 1:dim6
        bmchan = gammatoneFast(audio(:,ch,1,d4,d5,d6),cfs,fs,align);
        bm(:,ch,:,d4,d5,d6) = permute(bmchan,[1,3,2]);
                end
            end
        end
    end

    
    if isstruct(IN)
        OUT = IN; 
        OUT.audio = bm;
        OUT.bandID = cfs;
        OUT.funcallback.name = 'Gamma_tone_Fast_Hummersone2011.m';
        OUT.funcallback.inarg = {fs,hif,lof,bands,align};
    else
        OUT = bm;
    end
    varargout{1} = cfs;

else
    
    OUT = [];
end

%**************************************************************************
% Copyright (c) <YEAR>, <OWNER>
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
%  * Neither the name of the <ORGANISATION> nor the names of its contributors
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