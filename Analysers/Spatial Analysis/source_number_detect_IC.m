function OUT = source_number_detect_IC(IN,criterioncode)
% This function ports Choqueuse Vincent's 2009 function for the estimation
% of the number of sources in a multi-channel recording.
%
% Detection of the number of sources with information criteria. The method
% assumes that the noise is spatially white.
% Note that the number of sources cannot exceed the number of channels
%
% Reference: [WAX85] Wax, M. and Kailath, T., "Detection of signals by information
% theoretic criteria", IEEE Transactions on Acoustics, Speech and Signal
% Processing, 1985.

if nargin < 2 

    param = inputdlg({'Criterion: Akaike (0) or MDL (1)'},...
                      'Settings',... 
                      [1 60],... 
                      {'0'}); 

    param = str2num(char(param)); 

    if length(param) < 1, param = []; end 
    if ~isempty(param) 
        criterioncode = param(1);
    end
else
    param = [];
end
if isstruct(IN) 
    IN = choose_from_higher_dimensions(IN,2,1); 
    audio = IN.audio; % Extract the audio data

elseif ~isempty(param) || nargin > 1
    audio = IN;
end

if ~isempty(audio) && ~isempty(criterioncode)
    
    if size(audio,3) > 1
        audio = sum(audio,3); % mixdown if multiband
    end
    
    if criterioncode == 0
        criterion = 'AIC';
    else
        criterion = 'MDL';
    end
    
    % call the main function
    nb_sources = source_number_detection_IC(audio',criterion,1);

    if isstruct(IN)
        %OUT = IN; % You can replicate the input structure for your output
        OUT.nb_sources = nb_sources;
        OUT.funcallback.name = 'source_number_detect_IC.m';
        OUT.funcallback.inarg = {criterioncode};
    else
        OUT = nb_sources;
    end

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