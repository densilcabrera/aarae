function out = CentreTime_varTlim(data, fs, startthresh, highestband, lowestband, bpo, doplot)
% This function calculates centre time (ts) from an impulse
% response, using a range of integration times (from 10 ms to the end of
% the impulse response).
%
% This function is suitable for multi-channel and multi-band input (if
% required). If the input is multi-band, then the function does not apply
% any further bandpass filtering.
%
% Code by Densil Cabrera and Grant Cuthbert
% version 1.0 (15 December 2013)

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
        'Output data in table(s) (0|1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'-20','1','4000','125','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(answer)
        startthresh = str2num(answer{1,1});
        bpo = str2num(answer{2,1});
        highestband  = str2num(answer{3,1});
        lowestband  = str2num(answer{4,1});
        doplot = str2num(answer{5,1});
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

if ~isempty(IR) && ~isempty(fs) && ~isempty(startthresh) && ~isempty(highestband) && ~isempty(lowestband) && ~isempty(bpo) && ~isempty(doplot)
    %--------------------------------------------------------------------------
    % TRUNCATION
    %--------------------------------------------------------------------------

    % Check the input data dimensions
    [len, chans, bands] = size(IR);
    if bands > 1
        if isfield(data,'bandID')
            flist = data.bandID;
        else
            flist = 1:bands;
        end
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


    %--------------------------------------------------------------------------
    % FILTER AND CALCULATE CENTRE TIME
    %--------------------------------------------------------------------------


    % determine the filter frequencies
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
    end




    % filter and square
    if bands ==1
        bands = length(flist);
        if bpo ==3
            IR2 = thirdoctbandfilter(IR,fs,flist).^2;
        else
            IR2 = octbandfilter(IR,fs,flist).^2;
        end
    else
        IR2 = IR.^2;
    end

    IR2xt = IR2 .* repmat(((0:len-1)./fs)', [1,chans, bands]);


    %--------------------------------------------------------------------------
    % CALCULATE CENTRE TIME
    %--------------------------------------------------------------------------

    % List of integration limits
    tlim = 0.01:0.005:len/fs;

    Ts = zeros(length(tlim),chans,bands);
    for tcount = 1:length(tlim)

        Ts(tcount,:,:) = 1000 .* sum(IR2xt(1:floor(tlim(tcount).*fs),:,:),1) ...
            ./ sum(IR2(1:floor(tlim(tcount).*fs),:,:),1);
    end


    %--------------------------------------------------------------------------
    % OUTPUT STRUCTURE
    %--------------------------------------------------------------------------

    out.Ts = Ts;
    out.CentreFrequencies = flist;
    out.TimeLimits = tlim;
    out.Tsmax = squeeze(Ts(length(tlim),:,:));
    out.funcallback.name = 'CentreTime_varTlim.m';
    out.funcallback.inarg = {fs,startthresh,highestband,lowestband,bpo,doplot};

    %--------------------------------------------------------------------------
    % FIGURES
    %--------------------------------------------------------------------------


    [r, c] = subplotpositions(bands, 0.5);
    if doplot == 1, out.tables = []; end
    for ch = 1:chans

        figure('Name',['Centre time vs end truncation, channel ', num2str(ch)]);

        for b = 1:bands
            subplot(r,c,b)
            plot(tlim, Ts(:,ch,b),'Color',[0.5,0,0])
            title([num2str(flist(b)), ' Hz, ts = ', ...
                num2str(round(Ts(length(tlim),ch,b))), ' ms'])
            xlim([0 max(tlim)])
            xlabel('End Truncation Time (s)')
            ylabel('Centre Time (ms)')
            hold on
            plot([0, tlim(end)], ...
                [Ts(length(tlim),ch,b), Ts(length(tlim),ch,b)], ...
                'Color',[0.8,0,0], 'LineStyle', ':')
        end

        if doplot == 1
            % Table
            fig1 = figure('Name',['Centre Time ', num2str(ch)]);
            table1 = uitable('Data',permute(Ts(:,ch,:),[1,3,2]),...
                'ColumnName',num2cell(flist),...
                'RowName',num2cell(tlim*1000));
            set(table1,'ColumnWidth',{60});

            [~,datatable] = disptables(fig1,table1,{['Chan' num2str(ch) ' - Centre time']});
            out.tables = [out.tables datatable];
        end
    end
    if isstruct(data)
        if ~ismatrix(Ts)
            doresultleaf(Ts,'Centre Time [ms]',{'End_truncation_time'},...
                         'End_truncation_time', tlim,            's',          true,...
                         'channels',            chanID,          'categorical', [],...
                         'bands',               num2cell(flist), 'Hz',          false,...
                         'name','Centre_time');
        else
            doresultleaf(Ts,'Centre Time [ms]',{'End_truncation_time'},...
                         'End_truncation_time', tlim,   's',          true,...
                         'channels',            chanID, 'categorical', [],...
                         'name','Centre_time');
        end
    end
else
    out = [];
end
% eof

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
