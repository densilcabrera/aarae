function OUT = directivitydemo_ASM2016(IN,Bandfilter,Autocrop,Starttime,Endtime,AngleRange)
% This function plots and calculates the directivity (in one plane)
% recorded from a semi-circle or circle of microphones around a sound
% source.
%
% To use this function, first set up a semi-circular or circular microphone
% array centered the acoustic centre (or other reference point) of the
% source. The micrphones must have an even angular distribution over either
% 180 deg or 360 deg. The microphone channels should be in order of
% ascending angle, starting with channel 1. Channel 1 should be on the
% primary axis of interest in the measurement.
% If the source is a loudspeaker, usually you will want to measure an
% impulse response from it to each microphone (e.g. using a swept sinusoid
% test signal).
% Microphones should be calibrated (e.g. using AARAE's System Calibration
% via AARAE's recording window).
%
% If you don't have an anechoic room, you may still be able to make good
% measurements. This function facilitates the time-windowing of the direct
% sound, if the measurement has a sufficent reflection-free time period.
% Sometimes the measurement may be effectively done on a clear floor
% (reflective plane).
%


% notes on possible improvements:
% include a nominal distance input, and channel lag detection (frequency
% dependant) for inverse square law compensation?
%
% allow user to load microphone angles (from text file?)
%
% maybe display waveform analysed, so user can see if cropping is
% reasonable (does it include reflections?; does it cut at high amplitude?)
% use soft cropping, giving user control - potentially include interactive
% cropping loop
%
% allow user to specify filter bands
%
% do circular harmonic analysis to model the results

if size(IN.audio,2) == 1
    OUT = [];
    warndlg('Unable to analyse. This audio has only 1 channel but at least 2 are needed for directivity analysis.', 'Directivity analysis', 'modal')
    return
end

if nargin ==1
    
    param = inputdlg({'Broadband [0], Octave band [1] or 1/3-octave band [3]?';... % These are the input box titles in the
        'Autocrop start: No [0]; As ensemble [1]; Each channel individually [2]?';...
        'Start time (s)';...
        'End time (s)';
        'Mic angle range: 180 deg [1]; 360 deg [2]'},...
        'Directivity demo',...
        [1 60],...
        {'1';'1';'0';'0.01';'1'}); % Default values
    
    param = str2num(char(param));
    
    if length(param) < 5, param = []; end
    if ~isempty(param)
        Bandfilter = param(1);
        Autocrop = param(2);
        Starttime = param(3);
        Endtime = param(4);
        AngleRange = param(5);
    else
        % get out of here if the user presses 'cancel'
        OUT = [];
        return
    end
else
    param = [];
end


% *************************************************************************
if isstruct(IN)
    IN = choose_from_higher_dimensions(IN,3,1);

    if isfield(IN,'cal')
        IN = cal_reset_aarae(IN,0); % apply gain to set cal offset to 0 dB
    end
    
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    
    if size(audio,3)> 1
        Bandfilter = -1;
        % bandID is a vector, usually listing the centre frequencies of the
        % bands
        if isfield(IN,'bandID') % Get the band ID if it exists
            bandID = IN.bandID;
        else
            % asssign ordinal band numbers if bandID does not exist (as an
            % example of how to deal with missing data)
            bandID = 1:size(audio,3);
        end
    end
    
    % The name of the audio input could be useful in generating figures
    % (for example, in the title of a figure).
    if isfield(IN,'name') % Get the AARAE name if it exists
        name = IN.name; % this is a string
    else
        name = ''; % empty string - can be concatenated without any problems
    end
else
    %input must be an aarae audio structure
    OUT = [];
    return
end

% *************************************************************************
% Check that the required data exists for analysis to run
if ~isempty(audio) && ~isempty(fs)
    
    if Autocrop == 1
        % autocrop using ensemble peak
        audio = autocropstart_aarae(audio,-20,2);
    elseif Autocrop == 2
        % individual channel autocrop
        audio = autocropstart_aarae(audio,-20,1);
    end
    
    
    % the audio's length, number of channels
    [len,chans] = size(audio);
    Startindex = round(Starttime * fs) + 1;
    Endindex = round(Endtime * fs) + 1;
    if Endindex > len, Endindex = len; end
    if Startindex >= Endindex, Startindex = 1; end
    audio = audio(Startindex:Endindex,:,:);
    audio = [zeros(round(fs/100),chans,size(audio,3));audio;zeros(round(fs/100),chans,size(audio,3))];
    
    if Bandfilter == 1
        [audio,bandID] = octbandfilter_viaFFT(audio,fs,...
            [125,250,500,1000,2000,4000]);
    elseif Bandfilter == 3
        [audio,bandID] = thirdoctbandfilter_viaFFT(audio,fs,...
            [100,125,160,200,250,315,400,500,630,800,1000,1250,1600,...
            2000,2500,3150,4000,5000]);
    elseif ~exist('bandID','var')
        bandID = 0;
    end
    
    
    bands = size(audio,3);
    
    MeanSq = mean(audio.^2);
    MeanSq = MeanSq ./ repmat(MeanSq(1,1,:),[1,chans,1]); % normalise to ch 1
    L = pow2db(MeanSq);
    Lmax = max(max(max(L)));
    Lmin = min(min(min(L)));
    
    if AngleRange == 1
        theta = linspace(0,180,chans);
        if chans == 2
            Q2d = MeanSq(1,1,:) ./ repmat(mean(MeanSq,2),[1,chans,1]);
        else
            % double-count angles other than 0 and 180 deg in circle
            Q2d = MeanSq(1,1,:) ./ ...
                mean([mean(MeanSq,2);mean(MeanSq(1,2:end-1,:),2)],1);
        end
    else
        theta = linspace(0,360,chans+1);
        theta = theta(1:end-1);
        Q2d = MeanSq(1,1,:) ./ mean(MeanSq,2);
    end
    
    M = HSVplotcolours1(1, bands,[0.7 0.7]);
    % basic plot
    figure('Name',[name, ' sound level re 0 deg'])
    if Bandfilter == 3
        nsubplots = ceil(bands/3);
        [r, c] = subplotpositions(nsubplots, 0.75);
        n = 0;
        for b = 1:bands
            if mod(b,3)== 1
                n=n+1;
                subplot(r,c,n)
            end
            if mod(b,3)== 1
                linecolour = [1,0,0];
            elseif mod(b,3)== 2
                linecolour = [0,0.7,0];
            else
                linecolour = [0,0,1];
            end
            plot(theta,L(1,:,b),'Color',linecolour,...
                'DisplayName', [num2str(bandID(b)) ' Hz'])
            hold on
            if mod(b,3)== 1
                ylim([Lmin Lmax])
                if AngleRange == 1
                    xlim([0 180])
                else
                    xlim([0 360])
                end
                xlabel('Angle (deg)')
                ylabel('Level (dB)')
                %             title([num2str(bandID(b)), ' Hz, Q = ', num2str(round((100*Q2d(b)))/100),...
                %                 ' (DI = ', num2str(round(10*10*log10(Q2d(b)))/10), ' dB)'])
                grid on
            end
            %legend('show','Location','EastOutside')
        end
        for n = 1:nsubplots
            subplot(r,c,n)
            legend('show','Location','EastOutside')
        end
        
    else
        [r, c] = subplotpositions(bands, 0.75);
        for b = 1:bands
            subplot(r,c,b)
            plot(theta,L(1,:,b),'Color',M(b,:))
            ylim([Lmin Lmax])
            if AngleRange == 1
                xlim([0 180])
            else
                xlim([0 360])
            end
            xlabel('Angle (deg)')
            ylabel('Level (dB)')
            title([num2str(bandID(b)), ' Hz, Q = ', num2str(round((100*Q2d(b)))/100),...
                ' (DI = ', num2str(round(10*10*log10(Q2d(b)))/10), ' dB)'])
            grid on
        end
    end
    
    % polar plot
    try % the function polarplot was introduced by Matlab in 2016
        polarfig = figure('Name',[name, ' sound level re 0 deg']);
        for b = 1:bands
            if AngleRange == 1
                polarplot(pi*[theta theta(2:end)+180]/180,...
                    [L(1,:,b) flip(L(1,1:end-1,b),2)],...
                    'Color',M(b,:),...
                    'DisplayName', [num2str(bandID(b)) ' Hz'])
            else
                polarplot(pi*[theta 360]/180,...
                    [L(1,:,b) L(1,1,b)],...
                    'Color',M(b,:),...
                    'DisplayName', [num2str(bandID(b)) ' Hz'])
            end
            hold on
        end
        rlim([10*floor(Lmin/10) Lmax])
        title('Polar plot in dB relative to 0 deg')
        legend('show','Location','EastOutside')
    catch
        delete(polarfig)
        figure('Name',[name, ' sound magnitude re 0 deg']);
        for b = 1:bands
            if AngleRange == 1
                p=polar(pi*[theta theta(2:end)+180]/180,...
                    db2mag([L(1,:,b) flip(L(1,1:end-1,b),2)]));
                set(p,'Color',M(b,:),...
                    'DisplayName', [num2str(bandID(b)) ' Hz'])
            else
                p=polar(pi*[theta 360]/180,...
                    db2mag([L(1,:,b) L(1,1,b)]));
                set(p,'Color',M(b,:),...
                    'DisplayName', [num2str(bandID(b)) ' Hz'])
            end
            hold on
        end
        title('Polar plot in magnitude relative to 0 deg')
        legend('show','Location','EastOutside')
    end
    
    % plot of Q and DI as a function of frequency
    if bands > 1
        figure('Name',[name, ' Q and DI (planar)'])
        subplot(2,1,1)
        semilogx(bandID,permute(Q2d,[1 3 2]),'Marker','o')
        xlabel('Frequency (Hz)')
        ylabel('Directivity factor')
        grid on
        subplot(2,1,2)
        semilogx(bandID,permute(10*log10(Q2d),[1 3 2]),'Marker','o')
        xlabel('Frequency (Hz)')
        ylabel('Directivity index (dB)')
        grid on
    end
    
    % tables
    fig1 = figure('Name',['Directivity ' name]); % use the name if it is not empty
    table1 = uitable('Data',[Autocrop Starttime Endtime],...
        'ColumnName',{'Autocrop','Start (s)','End (s)'},...
        'RowName',{'Settings'});
    
    table2 = uitable('Data',[permute(Q2d,[1 3 2]);...
        10*log10(permute(Q2d,[1 3 2]));...
        permute(L,[2,3,1])],...
        'ColumnName',cellstr(num2str(bandID(:))),...
        'RowName',{char({'Q';'DI (dB)';...
        [repmat('L (dB) at ',[chans,1]) num2str(round(theta(:)*10)/10) repmat(' deg',[chans,1])]})});

    [~,tables] = disptables(fig1,[table1 table2]);
    OUT.tables = tables;
    
    OUT.funcallback.name = 'directivitydemo_ASM2016.m';
    OUT.funcallback.inarg = {Bandfilter,Autocrop,Starttime,Endtime,AngleRange};
    
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