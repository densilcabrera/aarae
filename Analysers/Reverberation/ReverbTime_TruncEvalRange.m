function OUT = ReverbTime_TruncEvalRange(IN, fs, bpo, highestband, lowestband, truncstep, evalrangestep, startlevel, startthresh, maketables, doplot)
% This function calculates apparent reverberation time as a function of
% end truncation time and evaluation range, using an impulse
% response as input.
%
% It may be useful in examining the reliability of reverberation time
% measurements, or in deciding how to best calculate reverberation time.
%
% The function can be relatively slow, because potentially it calculates a
% large number of apparent reverberation times.
% It generates a separate plot for each channel and band - be careful not
% to send it more channels/bands than you are interested in!
%
% The optional tables that it generates are not useful in themselves, but
% can be copied to the clipboard for pasting into a spreadsheet for further
% analysis.
%
% Code by Densil Cabrera
% version 1.01 (9 February 2014)

% -------------------------------------------------------------------------
% INPUTS AND SETTINGS
%
% The function can accept a structure as its first  and only input
% argument, in which case a dialog box is used to set parameter values. The
% input structure must have the .audio and .fs fields, and can have .bandID
% and .chanID fields (for multiband and multi-channel audio).
%
% Alternatively, the first input argument can be audio (as a vector or
% matrix), in which case all of the other input arguments must be used.
%
% 
% fs - audio sampling rate in Hz
%
% Bands per octave (bpo) allows for octave band (1) and 1/3-octave band (3)
% filtering. If set to 0, then no filtering is done. If the input is
% multi-band, then no filtering is done.
%
% Highest centre frequency (highestband) is the centre frequency of the
% highest octave or 1/3 octave band if bpo = 1 or 3 respectively. It is
% only used if filtering is done. Frequencies are automatically adjusted to
% IEC centre frequencies (and nominal frequencies are used in output data).
%
% Lowest centre frequency (lowestband) is the centre frequency of the
% lowest octave or 1/3 octave band if bpo = 1 or 3 respectively.
%
% The end of the impulse response will be truncated from 100 ms to its full
% length end in equally spaced steps. The truncation step (truncstep) is 
% this step. The smaller this value, the more detailed (and slow) the 
% analysis will be. The default value is 50 ms, which is probably suitable
% for most analyses.
%
% The end of the evaluation range is expressed in decibels, relative to the
% start of the Schroeder curve (not relative to the start of the evaluation
% range, which is usually -5 dB). Reverberation time is calculated from
% evaluation range end levels from -10 dB to -60 dB, in equal steps. The
% 'evaluation range step' (evalrangestep) is the size of this step -
% recommended values are 1, 2, 5 or 10 dB because these will enable the
% calculation of T20 and T30.
%
% The start of the evaluation range (startlevel) is normally -5 dB. 
% However, in some circumstances it may be useful to set it to 0 dB, or to 
% something else. The function does not allow start levels less than -7 dB.
%
% The impulse response start detection threshold (startthresh) is usually
% -20 dB (as recommended in ISO3382-1). This is used to remove parts of the
% wave prior to the start of the impulse response. This operation is
% performed for each channel. If multiband audio is input, then it is done
% for each band (but not if the audio is filtered within this function).
%
% The function allows tables to be generated, but typically the tables are
% too large to be displayed. Clicking on the grey areas between the tables
% copies the content to clipboard, which allows the data to be pasted into
% a spreadsheet for further analysis. Tables are in order from lowest to
% highest frequency band.

if nargin < 11, doplot = 0; end
if nargin < 10, maketables = 0; end
if nargin < 9, startthresh = -20; end
if nargin < 8, startlevel = -5; end
if nargin < 7, evalrangestep = 1; end
if nargin < 6, truncstep = 50; end
if nargin < 5, lowestband = 125; end
if nargin < 4, highestband = 8000; end
if nargin < 3, 
    bpo = 1;
    param = inputdlg({'Bands per octave (0 | 1 | 3)';...
        'Highest centre frequency (Hz)';...
        'Lowest centre frequency (Hz)';...
        'End truncation step (ms)';...
        'Evaluation range step (dB)';...
        'Evaluation range start level (dB)';...
        'Impulse response start detection threshold (dB)'; ...
        'Generate tables (0 | 1)';...
        'Generate plots (0 | 1)'},...
        'Settings',...
        [1 30],...
        {'1';'8000';'125';'50';'1';'-5';'-20';'0';'0'});
    
    param = str2num(char(param));
    
    if length(param) < 9, param = []; end
    if ~isempty(param)
        bpo = param(1);
        highestband = param(2);
        lowestband = param(3);
        truncstep = param(4);
        evalrangestep = param(5);
        startlevel = param(6);
        startthresh = param(7);
        maketables = param(8);
        doplot = param(9);
    else
        OUT = [];
        return
    end
else
    param = [];
end
if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1); 
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    if isfield(IN,'bandID')
        flist = IN.bandID;
    end
    if isfield(IN, 'chanID')
        chanID = IN.chanID;
    else
        chanID = cellstr([repmat('Chan',size(IN.audio,2),1) num2str((1:size(IN.audio,2))')]);
    end
elseif ~isempty(param) || nargin > 1
    audio = IN;
    
end

if ~isempty(audio) && ~isempty(fs) && ~isempty(bpo) && ~isempty(highestband) && ~isempty(lowestband) && ~isempty(truncstep) && ~isempty(evalrangestep) && ~isempty(startlevel)
    [len,chans,bands] = size(audio);
    
    
    % -------------------------------------------------------------------------
    % Find start point of IR
    % -------------------------------------------------------------------------
    
    
    % Preallocate
    m = zeros(1, chans,bands); % maximum value of the IR
    startpoint = zeros(1, chans,bands); % the auto-detected start time of the IR
    
    for dim2 = 1:chans
        for dim3 = 1:bands
            m(1,dim2,dim3) = max(audio(:,dim2,dim3).^2); % maximum value of the IR
            startpoint(1,dim2,dim3) = find(audio(:,dim2,dim3).^2 >= m(1,dim2,dim3)./ ...
                (10^(abs(startthresh)/10)),1,'first'); % Define start point
            
            if startpoint(1,dim2,dim3) >1
                
                % zero the data before the startpoint
                audio(1:startpoint(1,dim2,dim3)-1,dim2,dim3) = 0;
                
                % rotate the zeros to the end (to keep a constant data length)
                audio(:,dim2,dim3) = circshift(audio(:,dim2,dim3),-(startpoint(1,dim2,dim3)-1));
                
            end % if startpoint
        end
    end % for dim2
    
    
    % -------------------------------------------------------------------------
    % Filter signal if required, square and time-reverse
    % -------------------------------------------------------------------------
    
    
    % Determine frequencies of filters
    if bands == 1 && bpo ~= 0
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
        donotfilter = false;
    else
        if ~exist('flist', 'var')
            flist = 1:bands;
        end
        donotfilter = true;
    end
    
    % filter, square and time-reverse
    if ~donotfilter
        if bpo == 3
            % 1/3-octave band filterbank
            % decay = thirdoctbandfilter(audio,fs,flist).^2;
            
            order = [36,24];
            phasemode = 0;
            decay = thirdoctbandfilter_viaFFT(audio,fs,flist,order,0,1000,0,phasemode).^2;
        else
            % octave band filterbank
           % decay = octbandfilter(audio,fs,flist).^2;
            
            order = [12,12];
            phasemode = 0;
            decay = octbandfilter_viaFFT(audio,fs,flist,order,0,1000,0,phasemode).^2;
        end
    else
        decay = audio.^2;
    end
    
    bands = size(decay,3);
    
    
    % -------------------------------------------------------------------------
    % Calculate reverberation time
    % -------------------------------------------------------------------------
    
    
    % list of end truncation indices
    IRend = round(fs/10:truncstep*fs/1000:len);
    
    % list of evaluation range end indices
    ERendL = -10:-abs(evalrangestep):-65;
    
    % check validity of start level
    if startlevel < (ERendL(1)+3)
        startlevel = ERendL(1)+3'
        disp(['Evaluation range start level has been adjusted to ', num2str(startlevel), ' dB'])
    end
    
    Tstart = zeros(chans, bands);
    T = zeros(length(IRend),length(ERendL),chans, bands);
    
    for trunccount = 1:length(IRend)
        decaycurve = 10*log10(flipdim(cumsum(flipdim(decay(1:IRend(trunccount),:,:),1)),1));
        decaycurve = decaycurve - repmat(decaycurve(1,:,:),[IRend(trunccount),1,1]);
        for ch = 1:chans
            for bnd = 1:bands
                Tstart(ch, bnd) = find(decaycurve(:,ch,bnd) <= -5, 1, 'first');
            end
        end
        
        for ERcount = 1:length(ERendL)
            for ch = 1:chans
                for bnd = 1:bands
                    Tend = find(decaycurve(:,ch,bnd) <= ERendL(ERcount), 1, 'first');
                    % linear regression
                    p = polyfit((Tstart(ch,bnd):Tend)', ...
                        decaycurve(Tstart(ch,bnd):Tend,ch,bnd),1)';
                    % apparent reverberation time
                    T(trunccount,ERcount,ch,bnd) = 60/(-5 - ERendL(ERcount)) * ...
                        ((p(2)+ERendL(ERcount))/p(1) ...
                        -(p(2)-5)/p(1))/fs;
                end
            end
        end
    end
    
    
    % -------------------------------------------------------------------------
    % Plot
    % -------------------------------------------------------------------------
    T(isnan(T)) = 0;
    if maketables == 1, OUT.tables = []; end
    for ch =  1:chans
        if doplot == 1
            for bnd = 1:bands
                if exist('chanID', 'var')
                    figure('Name',['Apparent Reverberation time, ', num2str(chanID{ch,1})])
                else
                    figure('Name',['Apparent Reverberation time, channel ', num2str(ch)])
                end
                surfl(ERendL,IRend/fs,T(:,:,ch,bnd))
                xlabel('End of evaluation range (dB)')
                ylabel('Impulse response end truncation (s)')
                zlabel('Reverberation time (s)')


                if bands > 1
                    title([num2str(flist(bnd)), ' Hz'])
                end
            end
        
            if (startlevel == -5)
                T20index = find(ERendL == -25, 1, 'first');
                T30index = find(ERendL == -35, 1, 'first');
                if ~isempty(T20index) && ~isempty(T30index)
                    if exist('chanID', 'var')
                        figure('Name',['T20 & T30, channel ', num2str(chanID{ch,1})])
                    else
                        figure('Name',['T20 & T30, channel ', num2str(ch)])
                    end
                    subplot(2,1,1)

                    if bands > 1
                        c = selectcolormap('rainbow256');

                    end
                    for bnd = 1:bands
                        if bands > 1
                            cc = c(round(1+ 255 * (bnd-1)/(bands-1)),:);
                            else
                        cc = [1,0,0];
                        end
                        plot(IRend/fs,T(:,T20index,ch,bnd),...
                            'Color',cc,...
                            'DisplayName', [num2str(flist(bnd)), ' Hz'])
                        hold on
                    end
                    if bands > 1
                        legend('show','Location','EastOutside');
                    end
                    ylabel('Reverberation time (s)')
                    xlabel('Impulse response truncation point (s)')
                    title('T20')
                    hold off

                    subplot(2,1,2)
                    for bnd = 1:bands
                        if bands > 1
                            cc = c(round(1+ 255 * (bnd-1)/(bands-1)),:);
                            else
                        cc = [1,0,0];
                        end
                        plot(IRend/fs,T(:,T30index,ch,bnd),...
                            'Color',cc,...
                            'DisplayName', [num2str(flist(bnd)), ' Hz'])
                        hold on
                    end
                    if bands>1
                        legend('show','Location','EastOutside');
                    end
                    ylabel('Reverberation time (s)')
                    xlabel('Impulse response truncation point (s)')
                    title('T30')
                    hold off
                end
            end
        end
        if maketables == 1
            if ~isempty(T20index) && ~isempty(T30index)
                if exist('chanID', 'var')
                    fig1 = figure('Name',['Apparent Reverberation Time, channel ', num2str(chanID{ch,1})]);
                else
                    fig1 = figure('Name',['Apparent Reverberation Time, channel ', num2str(ch)]);
                end
                tables = [];
                for bnd = 1:bands
                    table1 = uitable('Data',T(:,:,ch,bnd),...
                        'ColumnName',num2cell(ERendL),...
                        'RowName',num2cell(IRend/fs));
                    tables = cat(2,tables,table1);
                end
                tablenames = cellstr([repmat(['Chan' num2str(ch) ' - '],length(flist),1) num2str(flist') repmat(' Hz',length(flist),1)]);
                [~,tables] = disptables(fig1,tables,tablenames);
                OUT.tables = [OUT.tables tables];
            end
        end
    end
    
    if isstruct(IN)
        T = permute(T,[2,1,3,4]);
        doresultleaf(T,'Reverb time [s]',{'End_truncation','End_eval_range'},...
                     'End_eval_range', ERendL,          'dB',          true,...
                     'End_truncation', IRend/fs,        's',           true,...
                     'channels',       chanID,          'categorical', [],...
                     'bands',          num2cell(flist), 'Hz',          false,...
                     'name','Reverb_time_TEC');
    end
    
    OUT.T = T;
    OUT.funcallback.name = 'ReverbTime_TruncEvalRange.m';
    OUT.funcallback.inarg = {fs,bpo,highestband,lowestband,truncstep,evalrangestep,startlevel};
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