function OUT = soundlogger(IN, windowtime, windowhop, tau, weighting, cal, fs)
% This function calculates sound level statistics in a series of windows
% over the duration of an audio signal (similar to a sound logger).
%
% Values returned are Lmax, L5, L10, L50, L90, Lmin and Leq.
%
% The audio wave can be weighted (A or C) or not.
%
% Temporal integration is applied to the rectified wave prior to analysis -
% use a time constant of 0.125 s for 'fast', 1 s for 'slow', or 0 for no
% integration.
%
% Code by Densil Cabrera
% version 0 (30 May 2014)

if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1); 
    if isfield(IN,'bandID')
        bandID = IN.bandID;
    end
    if isfield(IN,'chanID')
        chanID = IN.chanID;
    end
    if isfield(IN,'cal')
        cal = IN.cal;
    else
        cal = zeros(1,size(IN.audio,2));
        %disp('This audio signal has not been calibrated.')
    end
end

if nargin ==1
    
    param = inputdlg({'Window duration (s)';...
        'Hop size between windows (s)';...
        'Integration time constant (s)';
        'Audio weighting [a|c|z]'},...
        'Sound Logger Settings',...
        [1 60],...
        {'15';'15';'0.125';'a'}); % default values
    
    if length(param) < 4, param = []; end
    if ~isempty(param)
        windowtime = str2num(char(param(1)));
        windowhop = str2num(char(param(2)));
        tau = str2num(char(param(3)));
        weighting = char(param(4));
    else
        OUT = [];
        return
    end
else
    param = [];
end
if isstruct(IN)
    audio = IN.audio;
    fs = IN.fs;
    
    
elseif ~isempty(param) || nargin > 1
    
    audio = IN;
    
end


if ~isempty(audio) && ~isempty(fs)
    [len, chans, bands] = size(audio);
    
    WindowLength = round(windowtime*fs);
    Offset = round(windowhop*fs); % hopefully no need for rounding!
    nwin = round((len-WindowLength)/Offset); % number of windows
    if nwin < 1
        warndlg('Audio signal is shorter than the window length - unable to process with soundlogger')
        OUT = [];
        return
    end
    
    % apply calibration
    if exist('cal','var')
        if ~isempty(cal)
            if length(cal) == chans
                audio = audio .* repmat(10.^(cal./20),[len,1,bands]);
            elseif length(cal) == 1
                audio = audio .* 10.^(cal./20);
            end
        end
    end
    
    % apply weighting
    if weighting == 'a' || weighting == 'A'
        
        audio = Aweight(audio,fs);
    elseif weighting == 'c' || weighting == 'C'
        audio = Cweight(audio,fs);
    else
        weighting = 'z'; % change to 'z' if it is not already 'z'
    end
    
    % apply temporal integration
    % square and apply temporal integration
    if tau > 0
        % apply temporal integration so that percentiles can be derived
        % FILTER DESIGN
        E = exp(-1/(tau*fs)); % exponential term
        b = 1 - E; % filter numerator (adjusts gain to compensate for denominator)
        a = [1, -E];% filter denominator
        
        % rectify, integrate and square
        %audio=filter(b,a,abs(audio)).^2;
        
        % integrate the squared wave
        audio=filter(b,a,audio.^2);
        
    else
        % no temporal integration
        audio = audio.^2;
    end
    
    % analyse each window
    [Leq,Lmax,L5,L10,L50,L90,Lmin] = deal(zeros(nwin,chans,bands));
    for n = 1:nwin
        start = round((n-1)*Offset + 1);
        finish = start + WindowLength - 1;
        truncdata = audio(start:finish,:,:);
        
        Leq(n,:,:) = 10*log10(mean(truncdata));
        L = 10*log10(truncdata);
        Lmax(n,:,:) = max(L);
        L5(n,:,:) = prctile(L,95);
        L10(n,:,:) = prctile(L,90);
        L50(n,:,:) = median(L);
        L90(n,:,:) = prctile(L,10);
        Lmin(n,:,:) = min(L);
    end
    
    
    if ~exist('bandID','var')
        bandID = 1:bands;
    end
    
    
    t = windowhop*((1:nwin)-1)';
    ymax = 10*ceil(max(max(max(Lmax)))/10);
    if nwin > 1
        ymin = 10*floor(min(min(min(L90(2:end,:,:))))/10);
    else
        ymin = 10*floor(min(min(min(L90)))/10);
    end
    for ch = 1: chans
        if exist('chanID','var')
            chanstring = char(chanID(ch));
        else
            chanstring = ['ch ',num2str(ch)];
        end
        figure('Name',['Sound Logger Time Data ',chanstring])
        [r,c] = subplotpositions(bands,0.8);
        for b = 1:bands
            subplot(r,c,b)
            plot(t,Lmax(:,ch,b),...
                'Color',[1,0,0],...
                'DisplayName','Lmax');
            hold on
            plot(t,L5(:,ch,b),...
                'Color',[0.3,0.5,0],...
                'DisplayName','L5');
            plot(t,L10(:,ch,b),...
                'Color',[0,0.8,0],...
                'DisplayName','L10');
            plot(t,L50(:,ch,b),...
                'Color',[0,0,1],...
                'DisplayName','L50');
            plot(t,L90(:,ch,b),...
                'Color',[0.8,0,0.8],...
                'DisplayName','L90');
            plot(t,Lmin(:,ch,b),...
                'Color',[0.7,0.7,1],...
                'DisplayName','Lmin');
            plot(t,Leq(:,ch,b),...
                'Color',[1,0,0],...
                'LineStyle',':', ...
                'DisplayName','Leq');
            
            ylim([ymin ymax])
            xlabel('Window start time (s)')
            ylabel('Sound level (dB)')
            if bands > 1
                title([num2str(bandID(b)),' Hz'])
            end
            if bands == 1
                legend('Show','Location','EastOutside')
            end
        end
    end
    
    
    figure('Name','Sound Logger Distributions ')
    for b = 1:bands
        subplot(r,c,b)
        distributionPlot([reshape(Lmin(:,:,b),[nwin*chans,1]),...
            reshape(L90(:,:,b),[nwin*chans,1]), ...
            reshape(L50(:,:,b),[nwin*chans,1]),...
            reshape(L10(:,:,b),[nwin*chans,1]),...
            reshape(L5(:,:,b),[nwin*chans,1]),...
            reshape(Lmax(:,:,b),[nwin*chans,1]), ...
            reshape(Leq(:,:,b),[nwin*chans,1])]);
        set(gca,'XTickLabel',{'Lmin';'L90';'L50';'L10';'L5';'Lmax';'Leq'})
        ylim([ymin ymax])
        xlabel('Parameter')
        ylabel('Sound level (dB)')
        if bands > 1
            title([num2str(bandID(b)),' Hz'])
        end
    end
    
    if isstruct(IN)
        Log = cat(4,Lmin,L90,L50,L10,L5,Lmax,Leq);
        Log = permute(Log,[1,4,2,3]);
        dist = {'Lmin';'L90';'L50';'L10';'L5';'Lmax';'Leq'};
        if bands > 1
            doresultleaf(Log,'Sound level [dB]',{'time'},...
                         'time',         t,      's',           true,...
                         'distribution', dist,   'categorical', [],...
                         'channels',     chanID, 'categorical', [],...
                         'bands',        bandID, 'Hz',          false,...
                         'name','Sound_log');
        else
            doresultleaf(Log,'Sound level [dB]',{'time'},...
                         'time',         t,      's',           true,...
                         'distribution', dist,   'categorical', [],...
                         'channels',     chanID, 'categorical', [],...
                         'name','Sound_log');
        end
    end
    
    fig1 = figure('Name',['Sound Logger Summary ',weighting,'-weighted']);
    data = [windowtime;windowhop;tau];
    colname = {'Value'};
    rowname = {'Duration of each window (s)';...
        'Hop time from window to window (s)';...
        'Integration time constant (s)'};
    table1 = uitable('Data', data,...
        'ColumnName', colname,...
        'RowName', rowname);
    
    colname = {'Power mean';'Mean';'Median';'Minimum';'Maximum';'Range';...
        '90th prctile';'10th prctile';'90-10 percent range'};
    rowname = {'Leq';'Lmax';'L5';'L10';'L50';'L90';'Lmin'};
    table2 = [];
    if bands == 1, table2name = 'Soundlogger summary'; else table2name = cell(length(bands),1); end
    
    for b = 1:bands
        if bands > 1
            rowname = {['Leq ' num2str(bandID(b)) ' Hz'];...
                ['Lmax ' num2str(bandID(b)) ' Hz'];...
                ['L5 ' num2str(bandID(b)) ' Hz'];...
                ['L10 ' num2str(bandID(b)) ' Hz'];...
                ['L50 ' num2str(bandID(b)) ' Hz'];...
                ['L90 ' num2str(bandID(b)) ' Hz'];...
                ['Lmin ' num2str(bandID(b)) ' Hz']};
            table2name{b,1} = ['Soundlogger summary - ' num2str(bandID(b)) 'Hz band'];
        end
        data = [10*log10(mean(reshape(10.^(Leq(:,:,b)./10),[nwin*chans,1]))),...
            mean(reshape(Leq(:,:,b),[nwin*chans,1])),...
            median(reshape(Leq(:,:,b),[nwin*chans,1])),...
            min(reshape(Leq(:,:,b),[nwin*chans,1])),...
            max(reshape(Leq(:,:,b),[nwin*chans,1])),...
            max(reshape(Leq(:,:,b),[nwin*chans,1]))-min(reshape(Leq(:,:,b),[nwin*chans,1])),...
            prctile(reshape(Leq(:,:,b),[nwin*chans,1]),90),...
            prctile(reshape(Leq(:,:,b),[nwin*chans,1]),10),...
            prctile(reshape(Leq(:,:,b),[nwin*chans,1]),90)-prctile(reshape(Leq(:,:,b),[nwin*chans,1]),10);...
            ...
            10*log10(mean(reshape(10.^(Lmax(:,:,b)./10),[nwin*chans,1]))),...
            mean(reshape(Lmax(:,:,b),[nwin*chans,1])),...
            median(reshape(Lmax(:,:,b),[nwin*chans,1])),...
            min(reshape(Lmax(:,:,b),[nwin*chans,1])),...
            max(reshape(Lmax(:,:,b),[nwin*chans,1])),...
            max(reshape(Lmax(:,:,b),[nwin*chans,1]))-min(reshape(Lmax(:,:,b),[nwin*chans,1])),...
            prctile(reshape(Lmax(:,:,b),[nwin*chans,1]),90),...
            prctile(reshape(Lmax(:,:,b),[nwin*chans,1]),10),...
            prctile(reshape(Lmax(:,:,b),[nwin*chans,1]),90)-prctile(reshape(Lmax(:,:,b),[nwin*chans,1]),10);...
            ...
            10*log10(mean(reshape(10.^(L5(:,:,b)./10),[nwin*chans,1]))),...
            mean(reshape(L5(:,:,b),[nwin*chans,1])),...
            median(reshape(L5(:,:,b),[nwin*chans,1])),...
            min(reshape(L5(:,:,b),[nwin*chans,1])),...
            max(reshape(L5(:,:,b),[nwin*chans,1])),...
            max(reshape(L5(:,:,b),[nwin*chans,1]))-min(reshape(L5(:,:,b),[nwin*chans,1])),...
            prctile(reshape(L5(:,:,b),[nwin*chans,1]),90),...
            prctile(reshape(L5(:,:,b),[nwin*chans,1]),10),...
            prctile(reshape(L5(:,:,b),[nwin*chans,1]),90)-prctile(reshape(L5(:,:,b),[nwin*chans,1]),10);...
            ...
            10*log10(mean(reshape(10.^(L10(:,:,b)./10),[nwin*chans,1]))),...
            mean(reshape(L10(:,:,b),[nwin*chans,1])),...
            median(reshape(L10(:,:,b),[nwin*chans,1])),...
            min(reshape(L10(:,:,b),[nwin*chans,1])),...
            max(reshape(L10(:,:,b),[nwin*chans,1])),...
            max(reshape(L10(:,:,b),[nwin*chans,1]))-min(reshape(L10(:,:,b),[nwin*chans,1])),...
            prctile(reshape(L10(:,:,b),[nwin*chans,1]),90),...
            prctile(reshape(L10(:,:,b),[nwin*chans,1]),10),...
            prctile(reshape(L10(:,:,b),[nwin*chans,1]),90)-prctile(reshape(L10(:,:,b),[nwin*chans,1]),10);...
            ...
            10*log10(mean(reshape(10.^(L50(:,:,b)./10),[nwin*chans,1]))),...
            mean(reshape(L50(:,:,b),[nwin*chans,1])),...
            median(reshape(L50(:,:,b),[nwin*chans,1])),...
            min(reshape(L50(:,:,b),[nwin*chans,1])),...
            max(reshape(L50(:,:,b),[nwin*chans,1])),...
            max(reshape(L50(:,:,b),[nwin*chans,1]))-min(reshape(L50(:,:,b),[nwin*chans,1])),...
            prctile(reshape(L50(:,:,b),[nwin*chans,1]),90),...
            prctile(reshape(L50(:,:,b),[nwin*chans,1]),10),...
            prctile(reshape(L50(:,:,b),[nwin*chans,1]),90)-prctile(reshape(L50(:,:,b),[nwin*chans,1]),10);...
            ...
            10*log10(mean(reshape(10.^(L90(:,:,b)./10),[nwin*chans,1]))),...
            mean(reshape(L90(:,:,b),[nwin*chans,1])),...
            median(reshape(L90(:,:,b),[nwin*chans,1])),...
            min(reshape(L90(:,:,b),[nwin*chans,1])),...
            max(reshape(L90(:,:,b),[nwin*chans,1])),...
            max(reshape(L90(:,:,b),[nwin*chans,1]))-min(reshape(L90(:,:,b),[nwin*chans,1])),...
            prctile(reshape(L90(:,:,b),[nwin*chans,1]),90),...
            prctile(reshape(L90(:,:,b),[nwin*chans,1]),10),...
            prctile(reshape(L90(:,:,b),[nwin*chans,1]),90)-prctile(reshape(L90(:,:,b),[nwin*chans,1]),10);...
            ...
            10*log10(mean(reshape(10.^(Lmin(:,:,b)./10),[nwin*chans,1]))),...
            mean(reshape(Lmin(:,:,b),[nwin*chans,1])),...
            median(reshape(Lmin(:,:,b),[nwin*chans,1])),...
            min(reshape(Lmin(:,:,b),[nwin*chans,1])),...
            max(reshape(Lmin(:,:,b),[nwin*chans,1])),...
            max(reshape(Lmin(:,:,b),[nwin*chans,1]))-min(reshape(Lmin(:,:,b),[nwin*chans,1])),...
            prctile(reshape(Lmin(:,:,b),[nwin*chans,1]),90),...
            prctile(reshape(Lmin(:,:,b),[nwin*chans,1]),10),...
            prctile(reshape(Lmin(:,:,b),[nwin*chans,1]),90)-prctile(reshape(Lmin(:,:,b),[nwin*chans,1]),10)];...
            table2 = [table2 uitable('Data', data,...
            'ColumnName', colname,...
            'RowName', rowname)];
    end
    
    [~,tables] = disptables(fig1,[table1 table2],cat(1,{'Analysis details'},table2name));
    OUT.tables = tables;
    
    %     table1 = uitable('Data',[duration maximum minimum],...
    %                 'ColumnName',{'Duration','Maximum','Minimum'},...
    %                 'RowName',{'Results'});
    %     disptables(fig1,table1);
    % You may want to display your tables as barplots in teh AARAE
    % environment, in order to do this simply use the output of the
    % disptables funtion as follows:
    %       [~,table] = disptables(fig1,table1);
    % And include your table in the output data structure
    %       OUT.tables = table;
    % If you have multiple tables to combine in the figure, you can
    % concatenate them:
    %
    %       disptables(fig1,[table1 table2 table3]);
    %
    % You may export these tables to be displayed as bar plots as if you
    % were doing it for a single table:
    %       [~,tables] = disptables(fig1,[table1 table2 table3]);
    %       OUT.tables = tables;
    % The disptables function will take care of allocating each table to a
    % different barplot, there is no need to generate more than one .tables
    % field to display all your tables.
    
    % You may also include figures to display your results as plots.
    
    % All figures created by your function are stored in the AARAE
    % environment under the results box. If your function outputs a
    % structure in OUT this saved under the 'Results' branch in AARAE and
    % it's treated as an audio signal if it has both .audio and .fs fields,
    % otherwise it's displayed as data.
    % You may want to include your plots as part of the data variable
    % generated by AARAE, in order to do this use the getplotdata function
    % as follows:
    %       OUT.lines.myplot = getplotdata;
    % Use this function for as many charts as you want to include as output
    % from your function. Remember to call the getplotdata function after
    % you have designed your chart. Currently this function only supports
    % barplots and lines. E.g.:
    %
    %       plot(t,audio)
    %       OUT.lines.thischart = getplotdata;
    
    % And once you have your result, you should set it up in an output form
    % that AARAE can understand.
    
    %OUT = IN; % You can replicate the input structure for your output
    %OUT.audio = audio; % And modify the fields you processed
    % Or simply output the fields you consider necessary after
    % processing the input audio data, AARAE will figure out what has
    % changed and complete the structure. But remember, it HAS TO BE a
    % structure if you're returning more than one field:
    %   OUT = audio;  if you just want to return the audio,
    %   or,
    %   OUT.audio = audio; if you want to return two fields.
    %   OUT.fs = fs;
    OUT.funcallback.name = 'soundlogger.m'; % Provide AARAE
    % with the name of your function
    OUT.funcallback.inarg = {windowtime, windowhop, tau, weighting, cal, fs}; % assign all of the
    % input parameters that could be used to call the function
    % without dialog box to the output field param (as a cell
    % array) in order to allow batch analysing.
    
    % The processed audio data will be automatically displayed in AARAE's main
    % window as long as your output contains audio stored either as a single
    % variable: OUT = audio;, or it's stored in a structure along with any other
    % parameters: OUT.audio = audio;
else
    % AARAE requires that in case that the user doesn't input enough
    % arguments to generate audio to output an empty variable.
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