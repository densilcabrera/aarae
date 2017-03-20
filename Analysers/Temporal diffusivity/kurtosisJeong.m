function OUT = kurtosisJeong(IN,keepD, plotIR, stept, winlen)
% This function assesses the kurtosis of impulse response measurements as a
% measure of diffuseness according to:
%     C.H. Jeong, Kurtosis of room impulse responses as a diffuseness
%     measure for reverberation chambers, The Journal of the Acoustical
%     Society of America, vol. 139, no. 5, pp. 2833-2841, 2016.
% 
% Input audio can be multi-dimensional, e.g. input channels in dim2 and
% output channels in dim5.
%
% Users can change the default input parameters.
% The direct sound can be included in the calculation and the step time and
% windows size can be varied.
% 'Generate subplots' displays the normalised wave form and kurtosis of the
% user input window size.
% 

if nargin == 1
    param = inputdlg({'Keep Direct sound. No (0) is default, Yes (1) is optional';... % These are the input box titles in the
        'Generate subplots of IR and kurtosis (1 = yes, 0 = no)';...
        'Step time in s (default is 0.005 s)';...
        'window length in s (default is 0.02 s)'},...% inputdlg window.
        'User Inputs',... % This is the dialog window title.
        [1 30],... % You can define the number of rows per
        ...        % input box and the number of character
        ...        % spaces that each box can display at once
        ...        % per row.
        {'0';'1';'0.005'; '0.02'}); % And the preset answers for your dialog.
    
    param = str2num(char(param)); % Since inputs are usually numbers it's a
    % good idea to turn strings into numbers.
    % Note that str2double does not work
    % here.
    
    if length(param) < 4, param = []; end % You should check that the user
    % has input all the required
    % fields.
    if ~isempty(param) % If they have, you can then assign the dialog's
        % inputs to your function's input parameters.
        keepD = param(1);
        plotIR = param(2);
        stept = param(3);
        winlen = param(4);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
end

fs = IN.fs;

% frequency limits (125 Hz - 4 kHz octave band limits)
flo =88; %lower bound of 125 Hz Octave band
fhi = 5680; % upper bound of 4000 Hz octave band
order = [24,24];
phase = 0;
out = bandpass(IN, flo, fhi, order, fs, phase);
audio = out.audio;



% find combinations of output and input channels.
if size(IN.audio, 5) >1 && size(IN.audio,2) == 1
    audio = permute25_aarae(out);
    audio = audio.audio;
    if isfield(IN, 'OutchanID')
        chanID = IN.OutchanID'
    else
        chanID = []
    end
    
elseif size(IN.audio, 5) >1 && size(IN.audio,2) > 1
    % reshape dimensions
    ind = 1;
    filtstr = '[Cc]han';
    replace = '';
    for i = 1:size(IN.audio,5)
        for j = 1:size(IN.audio,2)
            tempaudio(:,ind) = audio(:,j,:,:,i);
            if isfield(IN, 'chanID') && isfield(IN, 'OutchanID')
                chanID{ind} = strcat(regexprep(['In',IN.chanID{j},' ' IN.OutchanID{i}], filtstr,replace));
            else
                chanID{ind} = ['chan: ', num2str(ind)];
            end
            ind = ind+1;
        end
    end
    audio = tempaudio;
else
    chanID = IN.chanID'
end

chanID2 = chanID;
chanID2{:,size(chanID,2)+1} = 'Av. all Ch';

% % % % % % commented out: kurtosis individual octave bands
% % % % % base = 10;
% % % % % test = 0;
% % % % % minfftlenfactor = 1000;
% % % % % zeropad =0;
% % % % %
% % % % % ob125 = octbandfilter_viaFFT(IN,fs,125,order,zeropad,minfftlenfactor,test,phase,base);
% % % % % ob250 = octbandfilter_viaFFT(IN,fs,250,order,zeropad,minfftlenfactor,test,phase,base);
% % % % % ob500 = octbandfilter_viaFFT(IN,fs,500,order,zeropad,minfftlenfactor,test,phase,base);
% % % % % ob1000 = octbandfilter_viaFFT(IN,fs,1000,order,zeropad,minfftlenfactor,test,phase,base);
% % % % % ob2000 = octbandfilter_viaFFT(IN,fs,2000,order,zeropad,minfftlenfactor,test,phase,base);
% % % % % ob4000 = octbandfilter_viaFFT(IN,fs,4000,order,zeropad,minfftlenfactor,test,phase,base);
% % % % % ob125 = ob125.audio;
% % % % % ob250 = ob250.audio;
% % % % % ob500 = ob500.audio;
% % % % % ob1000 = ob1000.audio;
% % % % % ob2000 = ob2000.audio;
% % % % % ob4000 = ob4000.audio;



[len,chans,bands] = size(audio);
endt = (.1+winlen)*fs;
if keepD ==1
    startt = 0;
else
    startt = 0.01*fs;
end
for i = 1:chans
    try
        peaks(i) = find(audio(:,i) == max(abs(audio(:,i))),1,'first');
    catch
        peaks(i) = find(audio(:,i) == -1* max(abs(audio(:,i))),1,'first');
    end
    starts(i) = find(abs(audio(:,i)) >= 0.1*abs(audio(peaks(i),i)),1,'first');
    starts(i) = starts(i)+startt;
    sigs(:,i) = audio(starts(i):starts(i)+endt,i);
    %normalise for plotting later
    sigs(:,i) = sigs(:,i)./max(abs(sigs(:,i)));
    % % % % %     sigsob125(:,i) = ob125(starts(i):starts(i)+endt,i);
    % % % % %     sigsob250(:,i) = ob250(starts(i):starts(i)+endt,i);
    % % % % %     sigsob500(:,i) = ob500(starts(i):starts(i)+endt,i);
    % % % % %     sigsob1000(:,i) = ob1000(starts(i):starts(i)+endt,i);
    % % % % %     sigsob2000(:,i) = ob2000(starts(i):starts(i)+endt,i);
    % % % % %     sigsob4000(:,i) = ob4000(starts(i):starts(i)+endt,i);
end

if length(sigs) < (.1+winlen)*fs
    timstr = num2str((.1+winlen)*1000);
    warndlg(['IR must contain ',timestr,' ms after peak'])
    out = [];
    return
else
    
    chunksize = winlen*fs;
    steptfs = stept*fs;
    
    % nchunks = 19; % for 100-120ms last chunk. nchunks = 15 for 80ms=100ms last chunk
    nchunks = (.1-startt/fs)./stept +1;
    for c = 1:chans
        for j = 1:nchunks
            chunkstart = (j-1)*steptfs +1;
            chunkend = chunkstart + chunksize-1;
            chunk = sigs(chunkstart:chunkend,c); 
            %         chunk = chunk./max(abs(chunk)); % normalise unesccessary as kurtosis function does it anyway
            % % % % %         chunk125 = sigsob125(chunkstart:chunkend,c);
            % % % % %         chunk250 = sigsob250(chunkstart:chunkend,c);
            % % % % %         chunk500 = sigsob500(chunkstart:chunkend,c);
            % % % % %         chunk1000 = sigsob1000(chunkstart:chunkend,c);
            % % % % %         chunk2000 = sigsob2000(chunkstart:chunkend,c);
            % % % % %         chunk4000 = sigsob4000(chunkstart:chunkend,c);
            jkurt(j,c) = kurtosis(chunk)-3;
            % % % % %         jkurt125(j,c) = kurtosis(chunk125)-3;
            % % % % %         jkurt250(j,c) = kurtosis(chunk250)-3;
            % % % % %         jkurt500(j,c) = kurtosis(chunk500)-3;
            % % % % %         jkurt1000(j,c) = kurtosis(chunk1000)-3;
            % % % % %         jkurt2000(j,c) = kurtosis(chunk2000)-3;
            % % % % %         jkurt4000(j,c) = kurtosis(chunk4000)-3;
        end
    end
    
    t = 0.01:stept:0.1;
    t = t(1:nchunks);
    if isfield(IN, 'name')
        fname = ['kurtosisJeong: ', IN.name];
    else
        fname = 'kurtosisJeong';
    end
    if plotIR == 0
        figure('Name' ,fname)
        plot(t,jkurt)
        xlim([0.01, .1])
        xlabel('Analysis window start (s)')
        ylabel(['kurtosis (',num2str(winlen*1000),' ms)'])
        legend(chanID, 'Location', 'southeastoutside')
    else
        figure('Name' ,fname)
        plot(t,jkurt)
        xlim([0.01, .1])
        xlabel('Analysis window start (s)')
        ylabel(['kurtosis (',num2str(winlen*1000),' ms)'])
        legend(chanID, 'Location', 'southeastoutside')
        % Generate plots for all source receievers
        % colums
        ncols = ceil(sqrt(size(sigs,2)));
        nrows = ceil(size(sigs,2)/ncols);
        if ncols*nrows == size(sigs,2) && size(sigs,2) >1
            ncols = ncols+1;
        end
        t2 = 0:1/fs:.12;
        figure('Name', fname, 'Position', [100,100,1250,1000])
        for n = 1:size(sigs,2)
            subplot(nrows, ncols, n)
            p1 =  plot(t2',sigs(:,n));
            hold on
            p2 = plot(t, jkurt(:,n));
            hold off
            xlim([0, .1+winlen])
            title(chanID{n})
            xlabel('time (s)')
        end  
        if n < ncols*nrows && n>1
            h =  subplot(nrows,ncols,n+1);
            title('Legend')
            p=get(h,'position');
            lh=legend(h,[p1;p2],'normalised waveform', ['kurtosis (', num2str(winlen*1000), 'ms)']);
            set(lh,'position',p, 'Units', 'normalized');
            axis(h,'off');
        else 
            legend('normalised waveform', ['kurtosis (', num2str(winlen*1000), 'ms)'], 'Location', 'southeast')
        end
    end
    
    meanjkurt = mean(jkurt,1);
    meanjkurtline = repmat(mean(meanjkurt), length(meanjkurt),1);
    
    figure('Name', fname)
    bar(meanjkurt, 'r')
    hold on
    bar(chans+1, mean(meanjkurt), 'b')
    ax = gca;
    ax.XTickLabelRotation =  90;
    ax.XTick = 1:length(chanID2);
    ax.XTickLabel = chanID2;
    plot(meanjkurtline,'--')
    hold off
    xlabel('Source Reveiver Combination');
    ylabel(['mean kurtosis from ', num2str(nchunks),' windows'])
    title(['mean kurtosis of ' num2str(size(meanjkurt,2)),' source & receiver combos = ',num2str(mean(meanjkurt))])
    
    
    % % % % %
    % % % % % meanjkurt125 = mean(jkurt125,1);
    % % % % % meanjkurt250 = mean(jkurt250,1);
    % % % % % meanjkurt500 = mean(jkurt500,1);
    % % % % % meanjkurt1000 = mean(jkurt1000,1);
    % % % % % meanjkurt2000 = mean(jkurt2000,1);
    % % % % % meanjkurt4000 = mean(jkurt4000,1);
    % % % % % M125 = [meanjkurt125 mean(meanjkurt125)];
    % % % % % M250 = [meanjkurt250 mean(meanjkurt250)];
    % % % % % M500 = [meanjkurt500 mean(meanjkurt500)];
    % % % % % M1000 = [meanjkurt1000 mean(meanjkurt1000)];
    % % % % % M2000 = [meanjkurt2000 mean(meanjkurt2000)];
    % % % % % M4000 = [meanjkurt4000 mean(meanjkurt4000)];
    % % % % %
    % % % % % f = figure('Name','octband kurtosis', ...
    % % % % %             'Position',[200 200 620 360]);
    % % % % %         dat1 = [M125;M250;M500;M1000;M2000;M4000;[meanjkurt, meanjkurtall(1)]];
    % % % % %         cnames1 = {IN.chanID{1:end,1},'average'}'; % fix this line. check chan ID's before continuing.
    % % % % %         rnames1 = {'125 Hz',...
    % % % % %             '250 Hz','500 Hz','1000 Hz','2000 Hz', '4000 Hz', 'Jeong (125 - 4000 Hz'};
    % % % % %         t1 =uitable('Data',dat1,'ColumnName',cnames1,'RowName',rnames1);
    % % % % %         set(t1,'ColumnWidth',{60});
    % % % % %         [~,tables] = disptables(f,[t1],{'KURTOSIS'});
    % % % % %         figure(4)
    % % % % %         bar(dat1)
    if isstruct(IN)
        OUT.funcallback.name = 'kurtosisJeong.m'; % Provide AARAE
        
        OUT.funcallback.inarg = {keepD,plotIR, stept, winlen};
    else
        out = [];
    end
    
end




%**************************************************************************
% Copyright (c) <2016, 2017>, <Jonothan Holmes>
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
%  * Neither the name of the <ORGANISATION> nor the names of its contributors
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
