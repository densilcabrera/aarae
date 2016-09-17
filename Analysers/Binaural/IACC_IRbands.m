function OUT = IACC_IRbands(IN,fs,startthresh,bpo,highestband,lowestband,maxlagtime,earlyendtime,doplot)
% This function calculates the interaural cross correlation function and
% IACC value from a binaural impulse response. Values are calculated in
% octave or 1/3-octave bands (or in the bands supplied by the audio if it
% has a 3rd dimension).
%
% The input audio should have two channels.
%
% If the input is multiband, then the bands-per-octave and highest and
% lowest centre frequency parameters are ignored.
%
% By default, the lag ranges between +- 1 ms, and this can be adjusted via
% the maximum lag time (maxlagtime).
%
% By default, the time at which early ends and late begins is 80 ms, and
% this can be adjusted (earlyendtime).
%
% If the function is run with octave band analysis, including the 500 Hz, 1
% kHz and 2 kHz band, then three-band average parameters are calculated,
% including the binaural quality index (BQI).
%
% When run outside of AARAE, all input arguments are required. Doplot can
% be use to suppress the output table and figure.
%
% Code by Densil Cabrera
% Version 1.02 (22 Feb 2016)

if nargin < 2
    param = inputdlg({'Threshold for IR start detection';...
        'Bands per octave (1 | 3)';...
        'Highest centre frequency (Hz)';...
        'Lowest centre frequency (Hz)';...
        'Maximum lag (ms)';...
        'Early/late time limit (ms)'},...
        'Settings',...
        [1 30],...
        {'-20';'1';'2000';'125';'1';'80'}); % Default values
    
    param = str2num(char(param));
    
    if length(param) < 6, param = []; end
    if ~isempty(param)
        startthresh = param(1);
        bpo = param(2);
        highestband = param(3);
        lowestband = param(4);
        maxlagtime = param(5);
        earlyendtime = param(6);
    end
else
    param = [];
end
if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1); 
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    if isfield(IN, 'bandID')
        flist = IN.bandID;
    end
    doplot = 1;
elseif ~isempty(param) || nargin > 1
    audio = IN;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

% check input dimensions
[len, chans, bands] = size(audio);

% zero-pad the audio if it is too short
if len < fs+1
    audio = [audio; zeros(fs+2-len,chans,bands)];
end

if chans > 2
    audio = audio(:,1:2,:);
    disp('Audio has more than 2 channels - only the first to channels will be analysed')
elseif chans == 1
    audio = repmat(audio,[1,2,1]);
    disp('This is a single channel audio wave - the wave has be replicated for a second channel')
end

if ~isempty(audio) && ~isempty(fs)
    
    % ---------------------------------------------------------------------
    % Find the start index and shift the audio to start at index 1
    % ---------------------------------------------------------------------
    
    % Preallocate
    m = zeros(1, chans, bands); % maximum value of the IR
    startpoint = zeros(1, chans, bands); % the auto-detected start index of the IR
    
    for dim2 = 1:chans
        for dim3 = 1:bands
            m(1,dim2) = max(audio(:,dim2,dim3).^2); % maximum value of the IR
            startpoint(1,dim2,dim3) = find(audio(:,dim2,dim3).^2 >= m(1,dim2,dim3)./ ...
                (10^(abs(startthresh)/10)),1,'first'); % Define start point
        end
    end
    
    startpoint = min(min(startpoint,[],3),[],2); % earliest startpoint
    if startpoint > 1
        audio(1:startpoint-1,:,:) = 0;
        audio = circshift(audio,-(startpoint-1));
    end % if startpoint
    
    
    
    % ---------------------------------------------------------------------
    % Filter the audio
    % ---------------------------------------------------------------------
    
    % Determine frequencies of filters
    if bands == 1
        hiband = round(10*log10(highestband));
        if hiband > 42, hiband = 42; end
        loband = round(10*log10(lowestband));
        if loband <15, loband = 15; end
        if hiband > loband
            if bpo == 3
                flist = 10.^((loband:hiband)./10);
            else
                flist = 10.^((loband:3:hiband)./10);
            end
        else
            flist = 10.^(hiband./10);
        end
        % convert exact frequencies to nominal
        flist = exact2nom_oct(flist);
    else
        if ~exist('flist', 'var')
            flist = 1:bands;
        end
    end
    
    
    % filter
    if bands == 1
        if bpo == 3
            % 1/3-octave band filterbank
            audio = thirdoctbandfilter(audio,fs,flist);
        else
            % octave band filterbank
            audio = octbandfilter(audio,fs,flist);
        end
    end
    
    
    % ---------------------------------------------------------------------
    % Cross-correlate
    % ---------------------------------------------------------------------
    
    maxlag = round(fs * maxlagtime / 1000);
    
    [cE,cL,cA] = deal(zeros(maxlag*2+1,bands));

    
    % Early
    Eend = floor(fs*earlyendtime/1000)+1;
    for bnd = 1:length(flist)
        [cE(:,bnd)] = xcorr(audio(1:Eend,1,bnd),audio(1:Eend,2,bnd),maxlag,'coeff');
    end
    ILDE = pow2db(mean(audio(1:Eend,1,:).^2) ./ mean(audio(1:Eend,2,:).^2));
    
    
    % Late
    Lstart = Eend+1;
    Lend = fs+1;
    for bnd = 1:length(flist)
        cL(:,bnd) = xcorr(audio(Lstart:Lend,1,bnd),audio(Lstart:Lend,2,bnd),maxlag,'coeff');
    end
    ILDL = pow2db(mean(audio(Lstart:Lend,1,:).^2) ./ mean(audio(Lstart:Lend,2,:).^2));
    
    % All
    for bnd = 1:length(flist)
        cA(:,bnd) = xcorr(audio(1:Lend,1,bnd),audio(1:Lend,2,bnd),maxlag,'coeff');
    end
    ILDA = pow2db(mean(audio(1:Lend,1,:).^2) ./ mean(audio(1:Lend,2,:).^2));
    
    IACC_E = max(abs(cE));
    IACC_L = max(abs(cL));
    IACC_A = max(abs(cA));
    
    % find lag times
    ind_E = zeros(length(flist),1);
    ind_L = zeros(length(flist),1);
    ind_A = zeros(length(flist),1);
    for bnd = 1:length(flist)
        ind_E(bnd) = find(abs(cE(:,bnd)) == max(abs(cE(:,bnd))), 1, 'first');
        ind_L(bnd) = find(abs(cL(:,bnd)) == max(abs(cL(:,bnd))), 1, 'first');
        ind_A(bnd) = find(abs(cA(:,bnd)) == max(abs(cA(:,bnd))), 1, 'first');
    end
    
    tau_E = (ind_E-maxlag) * 1000/fs;
    tau_L = (ind_L-maxlag) * 1000/fs;
    tau_A = (ind_A-maxlag) * 1000/fs;
    
    
    % ---------------------------------------------------------------------
    % Create output structure
    % ---------------------------------------------------------------------
    
    OUT.frequencies = flist;
    OUT.IACC_Early = IACC_E;
    OUT.IACC_Late = IACC_L;
    OUT.IACC_All = IACC_A;
    OUT.tau_Early = tau_E';
    OUT.tau_Late = tau_L';
    OUT.tau_All = tau_A';
    OUT.ILDE = permute(ILDE,[1,3,2]);
    OUT.ILDL = permute(ILDL,[1,3,2]);
    OUT.ILDA = permute(ILDA,[1,3,2]);
    
    % IACC E3 and L3
    ind500 = find(flist == 500, 1,'first');
    ind1000 = find(flist == 1000, 1,'first');
    ind2000 = find(flist == 2000, 1,'first');
    
    if ~isempty(ind500) && ~isempty(ind1000) && ~isempty(ind2000) && (ind2000 - ind500 == 2)
        OUT.IACC_E3 = mean(IACC_E(ind500:ind2000));
        OUT.IACC_L3 = mean(IACC_L(ind500:ind2000));
        OUT.BQI = 1-OUT.IACC_E3;
    end
    
    
    OUT.funcallback.name = 'IACC_IRbands.m';
    OUT.funcallback.inarg = {fs,startthresh,bpo,highestband,lowestband,maxlagtime,earlyendtime,doplot};
    
    % ---------------------------------------------------------------------
    % Display data in a table
    % ---------------------------------------------------------------------
    if doplot == 1
        OUT.tables = [];
        fig1 = figure('Name','Interaural Cross Correlation of a Room Impulse Response');
        table1 = uitable('Data',[OUT.IACC_Early; OUT.IACC_Late; OUT.IACC_All;...
            OUT.tau_Early;OUT.tau_Late;OUT.tau_All;...
            OUT.ILDE;OUT.ILDL;OUT.ILDA],...
            'ColumnName',num2cell(flist),...
            'RowName',{'IACC Early'; 'IACC Late'; 'IACC All';...
            'tau Early (ms)'; 'tau Late (ms)'; 'tau All (ms)';
            'ILD Early (dB)';'ILD Late (dB)';'ILD All (dB)'});
        set(table1,'ColumnWidth',{60});
        
        if isfield(OUT,'BQI')
            table2 = uitable('Data',[OUT.IACC_E3; OUT.IACC_L3; OUT.BQI],...
                'ColumnName',{'Result'},...
                'RowName',{'IACC E3'; 'IACC L3'; 'Binaural Quality Index'});
            
            [~,tables] = disptables(fig1,[table1 table2],{'IACC of an IR','IACC summary'});
            
        else
            [~,tables] = disptables(fig1,table1,{'IACC of an IR'});
        end
        OUT.tables = tables;
        % ---------------------------------------------------------------------
        % Display data in a figure
        % ---------------------------------------------------------------------
        
        figure('Name','Interaural Cross Correlation of a Room Impulse Response')
        [r, c] = subplotpositions(length(flist), 0.5);
        t = (-maxlag:maxlag) * 1000/fs;
        for bnd = 1:length(flist)
            subplot(r,c,bnd)
            plot(t',cE(:,bnd),'r');
            hold on
            plot([tau_E(bnd) tau_E(bnd)], [IACC_E(bnd) cE(ind_E(bnd),bnd)], 'ro:')
            
            plot(t',cL(:,bnd),'b');
            plot([tau_L(bnd) tau_L(bnd)], [IACC_L(bnd) cL(ind_L(bnd),bnd)], 'bo:')
            
            plot(t',cA(:,bnd),'k');
            plot([tau_A(bnd) tau_A(bnd)], [IACC_A(bnd) cA(ind_A(bnd),bnd)], 'ko:')
            
            text(-0.9,-0.5,['IACC(E) = ',num2str(round(100*IACC_E(bnd))/100)],'Color',[1,0,0]);
            text(-0.9,-0.7,['IACC(L) = ',num2str(round(100*IACC_L(bnd))/100)],'Color',[0,0,1]);
            text(-0.9,-0.9,['IACC(A) = ',num2str(round(100*IACC_A(bnd))/100)],'Color',[0,0,0]);
            
            xlim([t(1) t(end)])
            ylim([-1 1])
            
            xlabel('Lag (ms)')
            title([num2str(flist(bnd)), ' Hz'])
            hold off
        end
    end
    IACC = cat(3,cE,cL,cA);
    if isstruct(IN)
        doresultleaf(IACC,'Coefficient',{'Lag'},...
                     'Lag',       t,                      'ms',          true,...
                     'Frequency', num2cell(flist),        'Hz',          false,...
                     'IACC',      {'Early','Late','All'}, 'categorical', [],...
                     'name','IACC');
    end
else
    
    OUT = [];
end

%**************************************************************************
% Copyright (c) 2013-16, Densil Cabrera
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