function OUT = ChirpZplot(IN, f2, f1, r ,fs)
% This function implements Matlab's chirp z-transform and plots the absolute
% value of the result in decibels.
%
% r (radius) can be a single number (e.g., 1) or a row vector(e.g. 1:0.01:2)
%
% version 0 (22 March 2014)

if nargin ==1
    
    param = inputdlg({'Highest frequency (Hz)';...
        'Lowest frequency (Hz)';...
        'Radius (or radii), >=1'},...
        'Chirp z-transform',...
        [1 30],...
        {'20000';'100';'1'});
    
    
    if length(param) < 3, param = []; end
    if ~isempty(param)
        f2 = str2num(char(param(1)));
        f1 = str2num(char(param(2)));
        r = str2num(char(param(3)));
    end
else
    param = [];
end
if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1);
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    if isfield(IN,'bandID')
        bandID = IN.bandID;
    end
    if isfield(IN,'chanID')
        chanID = IN.chanID;
    end
    if isfield(IN,'cal')
        cal = IN.cal;
    else
        cal = zeros(1,size(audio,2));
        disp('This audio signal has not been calibrated.')
    end
elseif ~isempty(param) || nargin > 1
    audio = IN;
end


if ~isempty(audio) && ~isempty(fs)
    
    [len,chans,bands] = size(audio);
    if exist('cal','var')
        audio = audio .* repmat(10.^(cal./20),[len,1,bands]);
    end
    m=len;
    w = exp(-1i*2*pi*(f2-f1)/(m*fs));
    a = exp(1i*2*pi*f1/fs);
    
    
    S = zeros(m,chans,bands,length(r));
    for b = 1:bands
        for k=1:length(r)
            S(:,:,b,k) = czt(audio(:,:,b),m,w,r(k)*a);
        end
    end
    
    for ch = 1:chans
        if exist('chanID','var')
            chanstring = char(chanID(ch));
        else
            chanstring = ['chan ',num2str(ch)];
        end
        for b = 1:bands
            if exist('bandID','var')
                bandstring = num2str(bandID(b));
            else
                bandstring = num2str(b);
            end
            figure('Name',['Chirp z-transform ',chanstring,', ',bandstring])
            if size(S,4)>1
                imagesc(r,f1:(f2-f1)/(m-1):f2,10*log10(abs(squeeze(S(:,ch,b,:))).^2))
                set(gca,'YDir','normal');
                xlabel('Radius');
                ylabel('Frequency (Hz)')
                
            else
                plot(f1:(f2-f1)/(m-1):f2,10*log10(abs(S(:,ch,b)).^2))
                xlabel('Frequency (Hz)')
                ylabel('Level (dB)')
            end
        end
    end
    
    
    OUT.funcallback.name = 'ChirpZplot.m';
    OUT.funcallback.inarg = {f2, f1, r ,fs};
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