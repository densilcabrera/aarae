function out = ClarityandDefinition_varTlim(data, fs, startthresh, highestband, lowestband, bpo, maxtlim, doplot, endtime)
% This function calculates Clarity Index and Definition from an impulse
% response, using a range of integration times (from 10 ms to 100 ms, in 
% 5 ms steps).
%
% THIS IS QUITE A SLOW FUNCTION, because each evaluation range is filtered
% individually, as per ISO3382-1
%
% INPUT PARAMETERS:
% The threshold for IR start detection is recommended to be -20 dB in
% ISO3382-1. This is the level below the peak value, before the peak, at
% which the impulse response is assumed to start. This is detected on the
% unfiltered IR, and filtering is done after truncation.
%
% Either octave or 1/3-octave band filters can be applied. Multiband input
% is not accepted by this function because filtering must be done post
% truncation.
%
% The highest and lowest band frequencies can be entered. These are
% automatically adjusted to exact IEC band centre frequencies.
%
% In the case of an IR with a noisy tail, the end of the IR might have an
% effect on the results. Hence the end truncation point (in seconds) can be
% set. If the value is less than the duration of the wave, then the wave is
% truncated accordingly (otherwise the entire wave is used). The default
% value of 60 s is assumed to be longer than any likely input wave, and so
% results in no endpoint truncation. The endpoint time is the time from the
% start of the audio data (not the time from the detected IR start).
%
% Code by Densil Cabrera and Grant Cuthbert
% version 1.0 (15 December 2013)

if nargin < 8, endtime = 60; end
if nargin < 7, doplot = 1; end
if nargin < 6, bpo = 1; end
if nargin < 5, lowestband = 125; end
if nargin < 4, highestband = 4000; end
if nargin < 3
    startthresh = -20;
    %dialog box for settings
    prompt = {'Threshold for IR start detection', ...
        'Bands per octave (1 | 3)', ...
        'Highest band (Hz)', ...
        'Lowest band (Hz)', ...
        'End truncation time (s)', ...
        'Maximum early-late time (s)', ...
        'Plot (0|1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'-20','1','4000','125','60', '0.1','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        startthresh = str2num(answer{1,1});
        bpo = str2num(answer{2,1});
        highestband  = str2num(answer{3,1});
        lowestband  = str2num(answer{4,1});
        endtime  = str2num(answer{5,1});
        maxtlim = str2num(answer{6,1});
        doplot = str2num(answer{7,1});
    else
        out = [];
        return
    end
end
if isstruct(data)
    data = choose_from_higher_dimensions(data,3,1); 
    IR = data.audio;
    fs = data.fs;
else
    IR = data;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

if ~isempty(IR) && ~isempty(fs) && ~isempty(startthresh) && ~isempty(highestband) && ~isempty(lowestband) && ~isempty(bpo) && ~isempty(maxtlim) && ~isempty(doplot) && ~isempty(endtime)
    %--------------------------------------------------------------------------
    % TRUNCATION
    %--------------------------------------------------------------------------

    % Check the input data dimensions
    [len, chans, bands] = size(IR);
    if bands > 1
        warndlg('Input IR cannot be multiband','AARAE info','modal');
        out = [];
        return
    end
    if isfield(data,'chanID')
        chanID = data.chanID;
    end


    % Preallocate
    m = zeros(1, chans); % maximum value of the IR
    startpoint = zeros(1, chans); % the auto-detected start time of the IR

    for dim2 = 1:chans
        m(1,dim2) = max(IR(:,dim2).^2); % maximum value of the IR
        startpoint(1,dim2) = find(IR(:,dim2).^2 >= m(1,dim2)./ ...
            (10^(abs(startthresh)/10)),1,'first'); % Define start point

        if startpoint(1,dim2) > 1
            % zero the data before the startpoint
            IR(1:startpoint(1,dim2)-1,dim2) = 0;

            % rotate the zeros to the end (to keep a constant data length)
            IR(:,dim2) = circshift(IR(:,dim2),-(startpoint(1,dim2)-1));
        end
    end

    % End truncation
    if round(endtime*fs) < len
        IR = IR(1:round(endtime*fs),:);
    end


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


    %--------------------------------------------------------------------------
    % TRUNCATE, FILTER, SQUARE, SUM, AND CALCULATE CLARITY INDEX
    %--------------------------------------------------------------------------

    % List of temporal energy summation limits
    if maxtlim < 0.08, maxtlim = 0.08; end
    tlim = 0.01:0.005:maxtlim;
    index35 = find(tlim == 0.035, 1, 'first');
    index50 = find(tlim == 0.05, 1, 'first');
    index80 = find(tlim == 0.08, 1, 'first');

    % Pre-allocate
    ClarityIndex = zeros(length(tlim),chans,length(flist));

    % some zero-padding is applied to early, below, to avoid cutting off the
    % filters' decays
    for tcount = 1:length(tlim)
        if bpo == 3
            Early = sum(thirdoctbandfilter([IR(1:1+floor(fs*tlim(tcount)),:); ...
                zeros(round(fs/50),chans)],fs,flist).^2);
            Late = sum(thirdoctbandfilter(IR(ceil(fs*tlim(tcount)):end,:), ...
                fs,flist).^2);
        else
            Early = sum(octbandfilter([IR(1:1+floor(fs*tlim(tcount)),:); ...
                zeros(round(fs/50),chans)],fs,flist).^2);
            Late = sum(octbandfilter(IR(ceil(fs*tlim(tcount)):end,:), ...
                fs,flist).^2);
        end
        ClarityIndex(tcount,:,:) = 10*log10(Early./Late);
    end
    Definition = 10.^(ClarityIndex./10)./(10.^(ClarityIndex./10)+1);

    %--------------------------------------------------------------------------
    % OUTPUT STRUCTURE
    %--------------------------------------------------------------------------

    out.ClarityIndex = ClarityIndex;
    out.Definition = Definition;
    out.CentreFrequencies = flist;
    out.TimeLimits = tlim;
    out.C50 = ClarityIndex(index50,:,:);
    out.C80 = ClarityIndex(index80,:,:);
    out.funcallback.name = 'ClarityandDefinition_varTlim.m';
    out.funcallback.inarg = {fs,startthresh,highestband,lowestband,bpo,maxtlim,doplot,endtime};

    %--------------------------------------------------------------------------
    % FIGURES
    %--------------------------------------------------------------------------
    out.tables = [];
    if doplot == 1
        bands = length(flist);
        [r, c] = subplotpositions(bands, 0.5);

        for ch = 1:chans

            figure('Name',['Clarity Index, channel ', num2str(ch)]);

            for b = 1:bands
                subplot(r,c,b)
                plot(tlim*1000, ClarityIndex(:,ch,b),'Color',[0.5,0,0])
                title([num2str(flist(b)), ' Hz'])
                xlim([0 max(tlim)*1000])
                xlabel('Integration Limit (ms)')
                ylabel('Clarity Index (dB)')
                hold on
                plot(50,ClarityIndex(index50,ch,b),'ro')
                hold on
                plot([0 50],[ClarityIndex(index50,ch,b), ...
                    ClarityIndex(index50,ch,b)],':','Color', [1,0.8,0.8])
                plot(80,ClarityIndex(index80,ch,b),'ro')
                hold on
                plot([0 80],[ClarityIndex(index80,ch,b), ...
                    ClarityIndex(index80,ch,b)],':','Color', [1,0.8,0.8])
                text(10,ClarityIndex(index50,ch,b),...
                    [num2str(round(ClarityIndex(index50,ch,b)*10)/10), ' dB'])
                text(10,ClarityIndex(index80,ch,b),...
                    [num2str(round(ClarityIndex(index80,ch,b)*10)/10), ' dB'])
                hold off
            end

            % Table
            fig1 = figure('Name',['Clarity Index & Definition, channel ', num2str(ch)]);
            table1 = uitable('Data',permute(ClarityIndex(:,ch,:),[1,3,2]),...
                        'ColumnName',num2cell(flist),...
                        'RowName',num2cell(tlim*1000));
                    set(table1,'ColumnWidth',{60});

            table2 = uitable('Data',permute(Definition(:,ch,:),[1,3,2]),...
                'ColumnName',num2cell(flist),...
                'RowName',num2cell(tlim*1000));
            set(table2,'ColumnWidth',{60});

            [~,tables] = disptables(fig1,[table1 table2],{['Chan' num2str(ch) ' - Clarity index'],['Chan' num2str(ch) ' - Definition']});
            out.tables = [out.tables tables];
        end
        if isstruct(data)
            doresultleaf(ClarityIndex,'Clarity Index [dB]',{'Integration_limit'},...
                         'Integration_limit', tlim*1000,       'ms',          true,...
                         'channels',          chanID,          'categorical', [],...
                         'bands',             num2cell(flist), 'Hz',          false,...
                         'name','Clarity_index');
        end
        if chans == 1
            figure('Name', 'Clarity Index')
        else
            figure('Name', 'Clarity Index Average of All Channels')
        end

        meanClarity = 10*log10(mean(10.^(ClarityIndex./10),2));
        meanClarity = permute(meanClarity,[1,3,2]);
        All = [meanClarity(1:index35-1,:); ...
            meanClarity(index35+1:index50-1,:); ...
            meanClarity(index50+1:index80-1,:); ...
            meanClarity(index80+1:end,:)];

        semilogx(flist,All,'Color',[0.8 0.8 0.8])
        hold on
        semilogx(flist,meanClarity(index35,:),'-ro','LineWidth',1,'MarkerSize',8); 
        semilogx(flist,meanClarity(index50,:),'-r+','LineWidth',1,'MarkerSize',8);
        semilogx(flist,meanClarity(index80,:),'-rx','LineWidth',1,'MarkerSize',8);
        grid on
        xlabel('Band Frequency (Hz)');
        ylabel('Clarity Index (dB)');
        xlim([flist(1) flist(end)])
        title('Clarity index, 10-100 ms time limits (o C35; + C50; x C80)')
        %legend('show','Location','EastOutside');
        hold off
    end
else
    out = [];
end
end % eof

%**************************************************************************
% Copyright (c) 2013, Densil Cabrera and Grant Cuthbert
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
