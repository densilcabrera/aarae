function [OUT] = XFTmixingtime(IN, hoptime, maxtime, fs)
% This function visualizes the mixing time of a room impulse response
% through the eXtensive Fourier Transform, based on:
%
% G. Defrance and J.-D. Polack, 
% "Measuring the mixing time in auditoria,"
% Acoustics 2008, Paris.
%
% code by Densil Cabrera
% version 0 - beta (5 April 2014)

if nargin ==1

    param = inputdlg({'Hop size (s)';... 
                      'End duration (s)'},...
                      'Window title',... 
                      [1 30],... 
                      {'0.001';'0.5'}); 

    param = str2num(char(param)); 

    if length(param) < 2, param = []; end 
    if ~isempty(param) 
        hoptime = param(1);
        maxtime = param(2);
    end
else
    param = [];
end
if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1); 
    audio = IN.audio; 
    fs = IN.fs;       

    if isfield(IN,'chanID')
        chanID = IN.chanID;
    else
        chanID = cellstr([repmat('Chan',size(audio,2),1) num2str((1:size(audio,2))')]);
    end

    
    
    
elseif ~isempty(param) || nargin > 1
                       
    audio = IN;

end


if ~isempty(audio) && ~isempty(fs)
    [len, chans, bands] = size(audio);
    if bands > 1
        audio = sum(audio,3);
        disp('Multiband audio has been mixed prior to XFT analysis')
    end
    
    if maxtime+hoptime > len/fs
        maxtime = len/fs - hoptime;
    end
    
    nwin = round(maxtime/hoptime);
    
    audio = repmat(audio,[1,1,nwin]);
    
    hoplen = hoptime*fs;
    
    for n = 1:nwin
        audio(round(n*hoplen):end,:,n) = 0;
    end
    
    diplaylen = maxtime*fs;
    
    % unwraped phase in units of 2pi
    phase = unwrap(angle(fft(audio)))./(2*pi);
    
    D = zeros(nwin,chans);
    
    for ch = 1:chans
        for n = 1:nwin
            p = polyfit((1:len)', phase(:,ch,n),1)';
            % rms regression error
            D(n,ch) = (mean((((p(1)*(1:len)')+p(2))-phase(:,ch,n)).^2)).^0.5;
        end
    end
    
    
    
    %smoothed D
    runningavlen = 19;
    b = hann(runningavlen);
    b = b ./ sum(b);
    Dsmooth = filter(b,1,D);
    Dsmooth(1:round(runningavlen/2),:) = NaN;
    Dsmooth = circshift(Dsmooth,-round(runningavlen/2));
    
    
    
    for ch = 1:chans
        if exist('chanID','var')
            chanstring = char(chanID(ch));
        else
            chanstring = ['ch ',num2str(ch)];
        end
        figure('Name',['XFT results, ', chanstring])
        plot((0:diplaylen-1)'./fs,...
            audio(1:diplaylen,ch,nwin)./...
            max(abs(audio(1:diplaylen,ch,nwin))),'Color',[0.5,0.5,0.5])
        
        hold on
        plot((1:nwin)'*hoptime,...
            D(:,ch)./max(D(:,ch)),'r')
        
        plot((1:nwin)'*hoptime,...
            Dsmooth(:,ch)./max(D(:,ch)),'Color',[0.7,0,0.7])
        
        xlabel('Time (s)')
        ylabel('Normalized waveform and D')
        hold off
        %OUT.lines.(genvarname(['XFT_ch' num2str(ch)])) = getplotdata;
    end


    if isstruct(IN)
        All = cat(3,D./repmat(max(D),length(D),1),Dsmooth./repmat(max(D),length(D),1));
        doresultleaf(All,'normalized',{'Time'},...
                     'Time',     (1:nwin)'*hoptime,       's',           true,...
                     'channels', chanID,                  'categorical', [],...
                     'Mode',     {'Detailed','Smoothed'}, 'categorical', [],...
                     'name','Modal_density');
        %OUT = D; % This should be formatted to AARAE's generic data output
        % when the format has been defined.
        OUT.funcallback.name = 'XFTmixingtime.m';
        OUT.funcallback.inarg = {hoptime, maxtime};

    else
        
        OUT = [];
    end
%     varargout{1} = x;
%     varargout{2} = x;
%     varargout{3} = x;

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