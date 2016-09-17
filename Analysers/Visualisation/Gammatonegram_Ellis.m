function OUT = Gammatonegram_Ellis(IN,TWIN,THOP,N,FMAX,FMIN,USEFFT,WIDTH)
% This function calls Dan Ellis' gammatonegram, to produce a spectrogram
% inspired by auditory filters


if nargin ==1 

    param = inputdlg({'Window length (s)';...
        'Hop interval (s)';...
        'Number of filters';...
        'Maximum frequency (Hz)';...
        'Minimum frequency (Hz)';...
        'Use FFT (fast method) [0 | 1]';...
                      'Filter width (ERBs)'},...
                      'Settings',... % This is the dialog window title.
                      [1 30],... 
                      {'0.025';...
                      '0.01';...
                      '64';...
                      num2str(IN.fs/2);...
                      '50';...
                      '1';...
                      '1'}); 
    param = str2num(char(param)); 

    if length(param) < 7, param = []; end 
    if ~isempty(param) 
        TWIN = param(1);
        THOP = param(2);
        N = param(3);
        FMAX = param(4);
        FMIN = param(5);
        USEFFT = param(6);
        WIDTH = param(7);
    end
else
    param = [];
end
if isstruct(IN) 
    IN = choose_from_higher_dimensions(IN,2,1);
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    audio = sum(audio,3); % mixdown bands if present
    chans = size(audio,2);
elseif ~isempty(param) || nargin > 1
    disp('This function requires a structure as input')
return
end


if ~isempty(audio) && ~isempty(fs)
    if USEFFT == 1
            figure('Name','Gammatonegram - fast method')
        else
            figure('Name','Gammatonegram - accurate method')
    end
    
    % work out the size of the data for pre-allocation
    % using code from gammatonegram
    % (the efficiency of this can be greatly improved!)
    fcoefs = flipud(MakeERBFilters(fs, N, FMIN));
    XF = ERBFilterBank(audio(:,1),fcoefs);
    nwin = round(TWIN*fs);
    XE = XF.^2;
    hopsamps = round(THOP*fs);
    ncols = 1 + floor((size(XE,2)-nwin)/hopsamps);
    D = zeros(N,ncols,chans);
    for ch = 1:chans
        D(:,:,ch) = 20*log10(gammatonegram(audio(:,ch),fs,TWIN,THOP,N,FMIN,FMAX,USEFFT,WIDTH));
        
    end
    
    maxval = max(max(max(D)));
    
    for ch = 1:chans
        subplot(chans,1,ch)
        imagesc(D(:,:,ch)); axis xy
        caxis([maxval-60 maxval])
        colorbar
        if chans > 1
            title(['Channel ', num2str(ch)])
        end
    end
    OUT.funcallback.name = 'Gammatonegram_Ellis.m';
    OUT.funcallback.inarg = {TWIN,THOP,N,FMAX,FMIN,USEFFT,WIDTH};

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