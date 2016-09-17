function OUT = SpeechDistractionTester1(IN,method)
% This function is a testbed for prototype speech distraction index models
% Code by Densil Cabrera
% Version 0 (10 May 2016)

if nargin == 1
    param = inputdlg({'Method'},...
        'Speech Distraction Analysis',...
        [1 60],...
        {'1'}); % Default values
    
    param = str2num(char(param));
    
    if length(param) < 1, param = []; end
    if ~isempty(param)
        method = param(1);
    end
else
    OUT = [];
    return
end
if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1);
    audio = IN.audio;
    fs = IN.fs;
    if isfield(IN,'cal')
        IN = cal_reset_aarae(IN,0);
    end
    
    if isfield(IN, 'bandID')
        flist = IN.bandID;
    end
    if isfield(IN,'properties')
        if isfield(IN.properties,'name')
            inputname = IN.properties.name;
        else
            inputname = '';
        end
    else
        inputname = '';
    end
else
    OUT = [];
    disp('Unable to analyse: currently SpeechDistractionIndexTester1 is configured for structure input only')
    return
end

% check input dimensions
[len, chans, bands] = size(audio);

% zero-pad the audio if it is too short
if len < 2*fs+1
    audio = [audio; zeros(2*fs+2-len,chans,bands)];
    len = size(audio,1);
end

% apply calibration if available


% mixdown chans and bands for now (change later)
audio = mean(sum(audio,3),2);
chans = 1;
bands = 1;


if ~isempty(audio) && ~isempty(fs) && ~isempty(method)
    switch method
        otherwise
            % A-weighting
            audio = Aweight(audio,fs);
            
            % speech weighting (need to refine this)
            audio = LTASS_P50(audio,fs);
            
            % modulation analysis
            
            %Fm_all = [0.63,0.8,1,1.25,1.6,2,2.5,3.15,4,5,6.3,8,10,12.5];
            Fm_all = logspace(log10(0.1),log10(20),80);
            MTF = zeros(1,chans,bands,length(Fm_all));
            for j = 1:length(Fm_all)
                
                Intensity1 = audio.^2;
                Fm = Fm_all(j);
                % number of whole number cycles to use for each modulation frequency
                Fm_cycles = floor(len .* Fm./fs);
                
                % number of samples to use for each modulation frequency
                Fm_len = floor(fs.*Fm_cycles ./ Fm);
                
                % Time in seconds for each sample
                t=((1:len)-1)'./fs;
                
                % time replicated across channels and bands
                tch = repmat(t,[1,chans,bands]);
                
                
                MTF(1,:,:,j) = 2 * ((sum(Intensity1(1:Fm_len,:,:) ...
                    .* sin(2*pi*Fm*tch(1:Fm_len,:,:)))).^2 ...
                    + (sum(Intensity1(1:Fm_len,:,:) ...
                    .* cos(2*pi*Fm*tch(1:Fm_len,:,:)))).^2).^0.5 ...
                    ./ sum(Intensity1(1:Fm_len,:,:));
                
            end
            figure('name',inputname)
            semilogx(Fm_all',permute(MTF,[4,1,2,3]),'r')
            hold on
            smoothing = tukeywin(12,0.5);
            smoothing = smoothing ./ sum(smoothing);
            semilogx(Fm_all',filtfilt(smoothing,1,permute(MTF,[4,1,2,3])),'m')
            title(['mean = ' num2str(mean(mean(mean(mean(MTF)))))])
            xlabel('Modulation frequency (Hz)')
            ylabel('Modulation depth')
            xlim([Fm_all(1) Fm_all(end)])
            grid on
    end %switch method
    OUT.funcallback.name = 'SpeechDistractionTester1.m';
    OUT.funcallback.inarg = {method};
    
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