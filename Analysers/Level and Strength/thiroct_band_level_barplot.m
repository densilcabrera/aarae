function out = thiroct_band_level_barplot(in, fs, cal, showpercentiles, flo, fhi, tau, leveltype, dosubplots)
% This function generates 1/3-octave-band bar plots of various types of sound
% levels in decibels. There are four types of levels (determined by
% 'leveltype'):
%   * if leveltype == 0, Leq is plotted for each band, along with
%   percentiles
%   * if leveltype == 1, energy level is plotted for each band, along with
%   energy percentiles
%   * if leveltype == 2, Spectrum level Leq is plotted for each band, along
%   with percentiles. Spectrum level divides the power by the bandwidth of
%   each band (and so gives substantially lower values than the simple
%   octave band spectrum).
%   * if leveltype == 3, Spectrum energy level is plotted for each band,
%   along with percentiles.
%
% An integration time constant can be used (use 0.125 s for 'fast' and 1 s
% for 'slow'). If it has a value of 0, then temporal integration is not
% applied.
%
% Code by Densil Cabrera and Daniel Jimenez
% version 1.02 (19 October 2015)

if isstruct(in)
    %in = choose_from_higher_dimensions(in,3,1); 
    audio = in.audio;
    fs = in.fs;
    if isfield(in,'cal')
        cal = in.cal;
    else
        cal = 0;
        disp('This audio signal has not been calibrated.')
    end
    if isfield(in,'name') % Get the AARAE name if it exists
        name = in.name;
    else
        name = '';
    end
else
    audio = in;
    if nargin < 3
        cal = inputdlg({'Calibration offset (dB)'},...
                           'cal',1,{'0'});
        cal = str2num(char(cal));
    end
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
    name = '';
end
%if nargin < 8, dotables = 0; end
if nargin < 9, dosubplots = 0; end % default setting for multichannel plotting
if nargin < 8, leveltype = 0; end % power (use 1 for energy)
if nargin < 7, tau = 0; end % default temporal integration constant in seconds
if nargin < 6, fhi = 20000; end % default highest centre frequency
if nargin < 5, flo = 12.5; end % default lowest centre frequency
if nargin < 4
    showpercentiles = 1;
    param = inputdlg({'Individual plots [0] or subplots [1] for multiple channels';...
        'Integration time constant (s) - use 0.125 for fast, or 1 for slow, or 0 for no integration'; ...
        'Calibration offset (dB)'; ...
        'Highest third-octave band (Hz)'; ...
        'Lowest third-octave band (Hz)'; ...
        'Level type: band power [0], band energy [1], spectrum level power [2], spectrum level energy [3]';...
        'Show percentiles [0 | 1]'}, ...
        'Analysis and display parameters',1, ...
        {num2str(dosubplots);num2str(tau); num2str(cal);num2str(fhi); ...
        num2str(flo); num2str(leveltype); num2str(showpercentiles)});
    
   
            
    %if length(param) < 6, param = []; end
    if ~isempty(param)
        dosubplots = round(str2num(char(param(1))));
        tau = str2num(char(param(2)));
        cal = str2num(char(param(3)));
        fhi = str2num(char(param(4)));
        flo = str2num(char(param(5)));
        leveltype = str2num(char(param(6)));
        showpercentiles = str2num(char(param(7)));
    else
        out = [];
        return
    end
end

if ~isempty(audio) && ~isempty(fs)...
        && ~isempty(cal) && ~isempty(showpercentiles)...
        && ~isempty(flo) && ~isempty(fhi) && ~isempty(tau)...
        && ~isempty(dosubplots) && ~isempty(leveltype)


     [len,chans,bands,dim4,dim5,dim6]=size(audio);
    if bands > 1
        if isfield(in,'bandID')
            frequencies = in.bandID;
        else
            frequencies = 1:bands;
        end
        audiooct = audio;
    else
        % 1/3 octave band filterbank from AARAE: in Processors/Filterbanks
        % construct or reconstruct the structure
        hiband = round(10*log10(fhi));
        if hiband > 43, hiband = 42; end
        loband = round(10*log10(flo));
        if loband <11, loband = 11; end
        if hiband > loband
            flist = 10.^((loband:1:hiband)./10);
        else
            flist = 10.^(hiband./10);
        end
        % convert exact frequencies to nominal
        flist = exact2nom_oct(flist);
        [audiooct, frequencies] = thirdoctbandfilter_viaFFT(audio,fs,flist);
        frequencies = exact2nom_oct(frequencies);
        %disp(num2str(frequencies))
    end

    if tau > 0
        % apply temporal integration so that percentiles can be derived
        % FILTER DESIGN
        E = exp(-1/(tau*fs)); % exponential term
        b = 1 - E; % filter numerator (adjusts gain to compensate for denominator)
        a = [1 -E];% filter denominator
        % rectify, integrate and convert to decibels

        %Itemp=filter(b,a,abs(audiooct)).^2; % integrate the rectified wave
        Itemp=filter(b,a,audiooct.^2); % integrate the squared wave

    else
        % no temporal integration
        Itemp = audiooct.^2;
    end

    % apply calibration
    if length(cal) > 1
        for k = 1:chans
            Itemp(:,k,:,:,:,:) = 10.^((10*log10(Itemp(:,k,:,:,:,:)) + cal(k))./10);
        end
    else
        Itemp = 10.^((10*log10(Itemp) + cal)./10);
    end
    
    if leveltype == 1 || leveltype == 3
        Itemp = Itemp * len./fs;
        Lstring = 'Lenergy';
    else
        Lstring = 'Leq';
    end
    
    if leveltype == 2 || leveltype == 3
        for k = 1:length(frequencies)
            Itemp(:,:,k,:,:,:) = Itemp(:,:,k,:,:,:) ./ (frequencies(k).*2.^(0.5/3) - frequencies(k)./2.^(0.5/3));
        end
    end
    
    out.Leq = 10*log10(mean(Itemp));
    out.Leq = permute(out.Leq,[3,2,4,5,6,1]);
    
    out.Lmax = 10*log10(max(Itemp));
    out.Lmax = permute(out.Lmax,[3,2,4,5,6,1]);
    
    out.L1 = 10*log10(prctile(Itemp,99));
    out.L1 = permute(out.L1,[3,2,4,5,6,1]);
    
    out.L5 = 10*log10(prctile(Itemp,95));
    out.L5 = permute(out.L5,[3,2,4,5,6,1]);
    
    out.L10 = 10*log10(prctile(Itemp,90));
    out.L10 = permute(out.L10,[3,2,4,5,6,1]);
    
    out.L50 = 10*log10(median(Itemp));
    out.L50 = permute(out.L50,[3,2,4,5,6,1]);
    
    out.L90 = 10*log10(prctile(Itemp,10));
    out.L90 = permute(out.L90,[3,2,4,5,6,1]);
    
    out.funcallback.name = 'thiroct_band_level_barplot.m';
    out.funcallback.inarg = {fs,cal,showpercentiles,flo,fhi,tau,leveltype,dosubplots};
    
    ymax = 10*ceil(max(max(max(max(max(out.Lmax+5))/10))));
    ymin = 10*floor(min(min(min(min(min(out.L90))/10))));
    
    
    out.tables = [];
    for d4 = 1:dim4
        
        for d5 = 1:dim5
            
            for d6 = 1:dim6
                namestring = '';
                if dim4 > 1, namestring = ['Dim4 layer ' num2str(d4) ' ']; end
                if dim5 > 1, namestring = [namestring ', Dim5 layer ' num2str(d5) ' ']; end
                if dim6 > 1, namestring = [namestring ', Dim6 layer ' num2str(d6)]; end
                if dosubplots == 0
                    for ch = 1:chans
                        namestringch = [name ' Octave Band Spectrum, tau = ', ...
                            num2str(tau),' s, Channel ' num2str(ch) ', '...
                            namestring];
                        figure('name', namestringch)
                        
                        width = 0.5;
                        bar(1:length(frequencies),out.Leq(:,ch,d4,d5,d6),width,'FaceColor',[1,0.3,0.3],...
                            'EdgeColor',[0,0,0],'DisplayName', Lstring,'BaseValue',ymin);
                        hold on
                        
                        % x-axis
                        set(gca,'XTick',1:length(frequencies),'XTickLabel',num2cell(frequencies))
                        if (frequencies(1) ~= 1) && (frequencies(end) ~= length(frequencies))
                            xlabel('1/3-Octave Band Centre Frequency (Hz)')
                        else
                            xlabel('Band')
                        end
                        
                        % y-axis
                        if leveltype == 1 || leveltype == 3
                            ylabel('Energy Level (dB)')
                        else
                            ylabel('Level (dB)')
                        end
                        ylim([ymin ymax])
                        
                        if showpercentiles
                            plot(1:length(frequencies),out.Lmax(:,ch,d4,d5,d6),'Color',[0,0,0], ...
                                'Marker','o', 'DisplayName', 'Lmax')
                            hold on
                            plot(1:length(frequencies),out.L5(:,ch,d4,d5,d6),'Color',[0.1, 0.1, 0.1], ...
                                'Marker','o', 'LineStyle', ':', 'DisplayName', 'L5')
                            hold on
                            plot(1:length(frequencies),out.L10(:,ch,d4,d5,d6),'Color',[0.2, 0.2, 0.2], ...
                                'Marker','o', 'LineStyle', '--', 'DisplayName', 'L10')
                            hold on
                            plot(1:length(frequencies),out.L50(:,ch,d4,d5,d6),'Color',[0.3, 0.3, 0.3], ...
                                'Marker','o', 'LineStyle', '-.', 'DisplayName', 'L50')
                            hold on
                            plot(1:length(frequencies),out.L90(:,ch,d4,d5,d6),'Color',[0.4, 0.4, 0.4], ...
                                'Marker','o', 'LineStyle', ':', 'DisplayName', 'L90')
                            
                            legend('show','Location','EastOutside');
                            hold off
                        else
                            legend 'off'
                        end
                        
                        for k = 1:length(frequencies)
                            text(k-0.25,ymax-(ymax-ymin)*0.025, ...
                                num2str(round(out.Leq(k,ch,d4,d5,d6)*10)/10),'Color',[1,0.3,0.3])
                        end
                        fig = figure('Name',namestringch);
                        table1 = uitable('Data',[out.Leq(:,ch,d4,d5,d6),out.Lmax(:,ch,d4,d5,d6),out.L1(:,ch,d4,d5,d6),out.L5(:,ch,d4,d5,d6),out.L10(:,ch,d4,d5,d6),out.L50(:,ch,d4,d5,d6),out.L90(:,ch,d4,d5,d6)],...
                            'ColumnName',{Lstring,'Lmax','L1','L5','L10','L50','L90'},...
                            'RowName',num2cell(frequencies),'Parent',fig);
                        [~,tables] = disptables(fig,table1,{namestringch});
                        out.tables = [out.tables tables];
                    end
                else
                    [r, c] = subplotpositions(chans, 0.8);
                    figure('name', ['1/3-Octave Band Spectrum, tau = ', num2str(tau),' s'])
                    out.tables = [];
                    for ch = 1:chans
                        subplot(r,c,ch)
                        title(['Channel ', num2str(ch)])
                        
                        width = 0.5;
                        bar(1:length(frequencies),out.Leq(:,ch,d4,d5,d6),width,'FaceColor',[1,0.3,0.3],...
                            'EdgeColor',[0,0,0],'DisplayName', Lstring,'BaseValue',ymin);
                        hold on
                        
                        % x-axis
                        set(gca,'XTickLabel',num2cell(frequencies))
                        if (frequencies(1) ~= 1) && (frequencies(end) ~= length(frequencies))
                            xlabel('1/3-Octave Band Centre Frequency (Hz)')
                        else
                            xlabel('Band')
                        end
                        
                        % y-axis
                        if leveltype == 1 || leveltype == 3
                            ylabel('Energy Level (dB)')
                        else
                            ylabel('Level (dB)')
                        end
                        ylim([ymin ymax])
                        
                        if showpercentiles
                            plot(1:length(frequencies),out.Lmax(:,ch,d4,d5,d6),'Color',[0,0,0], ...
                                'Marker','o', 'DisplayName', 'Lmax')
                            hold on
                            plot(1:length(frequencies),out.L5(:,ch,d4,d5,d6),'Color',[0.1, 0.1, 0.1], ...
                                'Marker','o', 'LineStyle', ':', 'DisplayName', 'L5')
                            hold on
                            plot(1:length(frequencies),out.L10(:,ch,d4,d5,d6),'Color',[0.2, 0.2, 0.2], ...
                                'Marker','o', 'LineStyle', '--', 'DisplayName', 'L10')
                            hold on
                            plot(1:length(frequencies),out.L50(:,ch,d4,d5,d6),'Color',[0.3, 0.3, 0.3], ...
                                'Marker','o', 'LineStyle', '-.', 'DisplayName', 'L50')
                            hold on
                            plot(1:length(frequencies),out.L90(:,ch,d4,d5,d6),'Color',[0.4, 0.4, 0.4], ...
                                'Marker','o', 'LineStyle', ':', 'DisplayName', 'L90')
                            if ch == chans
                                legend('show');
                            else
                                legend 'off'
                            end
                            hold off
                        else
                            legend 'off'
                        end
                        
                       
                        
                    end
                    for ch = 1:chans
                        namestringch = [name ' 1/3-Octave Band Spectrum, tau = ', ...
                            num2str(tau),' s, Channel ' num2str(ch) ', '...
                            namestring];
                        fig = figure('Name',namestringch);
                        table1 = uitable('Data',[out.Leq(:,ch,d4,d5,d6),out.Lmax(:,ch,d4,d5,d6),out.L1(:,ch,d4,d5,d6),out.L5(:,ch,d4,d5,d6),out.L10(:,ch,d4,d5,d6),out.L50(:,ch,d4,d5,d6),out.L90(:,ch,d4,d5,d6)],...
                            'ColumnName',{Lstring,'Lmax','L1','L5','L10','L50','L90'},...
                            'RowName',num2cell(frequencies),'Parent',fig);
                        [~,tables] = disptables(fig,table1,{namestringch});
                        out.tables = [out.tables tables];
                    end
%                     tablenames = cellstr([repmat('Chan',chans,1),num2str((1:chans)'),repmat(' - Octave band spectrum',chans,1)]);
%                     [~,out.tables] = disptables(fig,out.tables,tablenames);
                end
            end
        end
    end
else
    out = [];
end

%**************************************************************************
% Copyright (c) 2013-2015, Densil Cabrera & Daniel Jimenez
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