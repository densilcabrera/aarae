function output=TN_LeeCabreraMartensJASA2012v3(RIR, level, ERdB_User)

% This function calculates the loudness-based reverberation time (TN) and
% early decay time (EDTN).
%
% EDTN and TN are similar to their conventional counterparts, but the
% loudness decay envelope is used instead of the time-reversed squared
% sound pressure envelope. Lee et al. (2010 and 2012) reported that the
% loudness-based parameters outperform the conventional EDT and RT
% for moderate listening levels and reverberation conditions.
%
% The loudness decay function of a room impulse response is calculated
% following Moore,Glasberg & Baer (1997)
% with Glasberg & Moore (2002) dynamic loudness
% and Moore & Glasberg (2007) binaural loudness
%
% FOR MORE INFORMATION, SEE:
% D. Lee, D. Cabrera and W.L. Martens (2012) "The effect of loudness
% on the reverberance of music: Reverberance prediction using loudness
% models," Journal of the Acoustical Society of America 131(2), 1194-1205.
%
% D. Lee and D.Cabrera (2010) "Effect of listening level and background
% noise on the subjective decay rate of room impulse responses: Using
% time-varying loudness to model reverberance," Applied Acoustics 71,
% 801-811.
%
% To run this code, the signal processing and curve fitting toolboxes
% are necessary.
%
% By Doheon Lee and Densil Cabrera (2013)
% version 3.00 (15 Sep 2016)

% INPUT
%     RIR:
%       -a room imuplse response (in a form of matrix or wav file) or a
%       structure including the fields .audio and .fs
%
%     LEVEL:
%       - EDTN and TN are sensitive to the playback level of signals, so
%           this input should be carefully determined.
%       - When users derive EDTN and TN for running stimuli (e.g., speech
%             or music), LEVEL should be a LAeq, at which the stimulus is
%             expected to be played in a room where the RIR is meausred. LAeq
%             is the A-weighted sound pressure level averaged over time.
%       - For an inpulsive signal, LEVEL should be a expected LAFmax level of
%           that impulsive signal. LAFmax is the maximum A-weighted sound
%           pressure level with a fast temporal integratoin of 125 ms.
%
%     ERdB:
%       - An evaluation range in dB. Default settings are [0 -10] for EDTN
%         and [-5 -25] for TN. It is recommended to use the default
%         setting, unless users have
%
% OUTPUT
%       - EDTN,TN
%       - TN and TN_User values
%       - Instantaneous Loudness, Short-term Loudness, Long-term Loudness

% References:
%  Lee, D., and Cabrera, D. (2010). Appl. Acoust. 71, 801-811.
%  Lee, D., Cabrera, D., and Martens, W. L. (2012). J. Acoust. Soc. Am. 131, 1194-1205.
%  Moore, B.J.C., Glasberg, B.R., and Baer, T. (1997). J. Acoust. Soc. Am. 45, 224-240.
%  Glasberg, B.R., and Moore, B.C.J. (2002). J. Audio. Eng. Soc. 50, 331-342.
%  Moore, B.J.C., and Glasberg, B.R. (2007). J. Acoust. Soc. Am. 121, 1604-1612

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2013, Doheon Lee and Densil Cabrera
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
%  * Neither the name of the University of Sydney nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Load wave data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(RIR)
    data = RIR.audio;
    fs = RIR.fs;
    
    if length(size(data))>2
        data = mean(data,3);
        disp('Multiband data has been mixed together.')
    end
    if nargin < 2
        prompt = {'Listening Level (dBA):', ...
            'TN_User Evaluation Range Start Point (~dB):', ...
            'TN_User Evaluation Range End Point (~dB):'};
        dlg_title = 'Settings';
        num_lines = 1;
        def = {'70','-5','-15'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if isempty(answer)
            output = [];
            return
        else
            level = str2num(answer{1,1});
            TN_User_start = str2num(answer{2,1});
            TN_User_end = str2num(answer{3,1});
        end
    end
    
elseif ischar(RIR)== 1
    [data, fs]=audioread(char(RIR));
else
    data=RIR;
    fs=str2double(input('Sampling Rate in Hz: ', 's'));
end

% USE FOR TESTING THE OUTPUT OF MGBLoudness FOR 1K AT 40 DECIBEL
% IMPORTANT: TESTING THIS WITH THE FOURTH INPUT VARIABLE = 'Leq';
% t=(1:fs);
% time=(t-1)./fs;
% data=sin(2*pi*1000*time);
% level=40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% EDTN AND TN CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LoudnessInstan, LoudnessShort, LoudnessLong,LoudnessTime]=MGBLoudness(data, fs, 2, level, 'LAeq', 1, 500, 0);

% SMOOTHING THE LOUDNESS DECAY CURVE AND DETECTING THE LOUDNESS OF DIRECT SOUND
[bb,a]=butter(3, 50/1000, 'low'); % note that the sampling rate of the output of 'MGBLoundssforTN' is 1000 Hz.
LoudnessDecayEnvelope=filtfilt(bb,a,LoudnessShort); clear bb a
pks=findpeaks(LoudnessDecayEnvelope, 'minpeakheight', max(LoudnessDecayEnvelope)*0.9);
LoudnessDirectSound=pks(1); % Peak loudness (presumably correspodning to the direct sound)
LoudnessDecaylog=log(LoudnessDecayEnvelope);

% % WITHOUT SMOOTHING
% LoudnessDecayEnvelope=LoudnessShort;
% LoudnessDirectSound=max(LoudnessDecayEnvelope);
% LoudnessDecaylog=log(LoudnessDecayEnvelope);

% EDTN CALCULATION
ERdB = [0 -10]; % Evaluation range in decibel-like units
ERsone=(10.^(ERdB/20)).^0.6;%Evaluation range in sone
Point1_EDTN= find(LoudnessDecayEnvelope==(LoudnessDirectSound*ERsone(1)),1,'last');
Point2_EDTN= find(LoudnessDecayEnvelope>=(LoudnessDirectSound*ERsone(2)),1,'last');
P_EDTN=polyfit(LoudnessTime(Point1_EDTN:Point2_EDTN),LoudnessDecaylog(Point1_EDTN:Point2_EDTN),1);
BestFitLineEDTN=polyval(P_EDTN,LoudnessTime);

if Point2_EDTN >= length(LoudnessDecayEnvelope)
    disp('Evaluation range is beyond the length of the impulse response')
    output.EDTN = [];
else
    output.EDTN=log((10^(-60/20))^0.6)./P_EDTN(1);
end

% TN CALCULATION
ERdB = [-5 -25]; % Evaluation range in decibel-like units
ERsone=(10.^(ERdB/20)).^0.6;%Evaluation range in sone
Point1_TN= find(LoudnessDecayEnvelope>=(LoudnessDirectSound*ERsone(1)),1,'last');
Point2_TN= find(LoudnessDecayEnvelope>=(LoudnessDirectSound*ERsone(2)),1,'last');
P_TN=polyfit(LoudnessTime(Point1_TN:Point2_TN),LoudnessDecaylog(Point1_TN:Point2_TN),1);
BestFitLineTN=polyval(P_TN,LoudnessTime);

if Point2_TN >= length(LoudnessDecayEnvelope)
    disp('The evaluation range is beyond the length of the impulse response')
    output.TN = [];
else
    output.TN=log((10^(-60/20))^0.6)./P_TN(1);
end

% CALCULATE TN_USER (USER-SPECIFIED EVALUATION RANGE)
if exist('TN_User_start', 'var') || nargin == 3
    if exist('TN_User_start', 'var')
        ERdB = [TN_User_start TN_User_end]; % Evaluation range in decibel-like units
    else
        ERdB = ERdB_User;
    end
    ERsone=(10.^(ERdB/20)).^0.6;%Evaluation range in sone
    if ERdB(1) == 0
        Point1_TNU= find(LoudnessDecayEnvelope==(LoudnessDirectSound*ERsone(1)),1,'last');
        Point2_TNU= find(LoudnessDecayEnvelope>=(LoudnessDirectSound*ERsone(2)),1,'last');
    else
        Point1_TNU=find(LoudnessDecayEnvelope>=(LoudnessDirectSound*ERsone(1)),1,'last');
        Point2_TNU= find(LoudnessDecayEnvelope>=(LoudnessDirectSound*ERsone(2)),1,'last');
    end
    
    P_TNU=polyfit(LoudnessTime(Point1_TNU:Point2_TNU),LoudnessDecaylog(Point1_TNU:Point2_TNU),1);
    BestFitLineTNU=polyval(P_TNU,LoudnessTime);
    
    if Point2_TNU >= length(LoudnessDecayEnvelope)
        disp('Evaluation range is beyond the length of the impulse response')
        output.TN_User = [];
    else
        output.TN_User=log((10^(-60/20))^0.6)./P_TNU(1);
    end
    output.funcallback.name = 'TN_LeeCabreraMartensJASA2012.m';
    output.funcallback.inarg = {level,ERdB};
end

output.LoudnessShort=LoudnessShort;
output.LoudnessLong=LoudnessLong;
output.LoudenssInstan=LoudnessInstan;
output.LoudnessTime=LoudnessTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('name','Loudness-based Reverberance Parameters')
if isfield(output, 'TN_User')
    subplot(1,3,1)
else
    subplot(1,2,1)
end

warning ('off', 'MATLAB:Axes:NegativeDataInLogAxis');
semilogy(LoudnessTime, LoudnessShort, LoudnessTime, LoudnessDecayEnvelope, ':k', ...
    [LoudnessTime(Point1_EDTN) LoudnessTime(Point2_EDTN)], ...
    exp([BestFitLineEDTN(Point1_EDTN) BestFitLineEDTN(Point2_EDTN)]), '--xr');
hold on
title(['EDTN: ', num2str(output.EDTN), ' s'])
xlabel('Time (s)', 'fontsize', 15);
ylabel('Loudness (sone)', 'fontsize', 15);
ymax = 10.^ceil(log10(max([max(LoudnessShort) max(BestFitLineEDTN)])));
ylim([ymax/1000 ymax])
%axis([0 3 1 4]);
if isfield(output, 'TN_User')
    subplot(1,3,2)
else
    subplot(1,2,2)
end

semilogy(LoudnessTime, LoudnessShort, LoudnessTime, LoudnessDecayEnvelope, ':k', ...
    [LoudnessTime(Point1_TN) LoudnessTime(Point2_TN)], ...
    exp([BestFitLineTN(Point1_TN) BestFitLineTN(Point2_TN)]), '--xr');
title(['TN: ', num2str(output.TN), ' s'])
xlabel('Time (s)', 'fontsize', 15);
ylabel('Loudness (sone)', 'fontsize', 15);
ymax = 10.^ceil(log10(max([max(LoudnessShort) max(BestFitLineTN)])));
ylim([ymax/1000 ymax])
%axis([0 3 1 4]);

if isfield(output, 'TN_User')
    subplot(1,3,3)
    semilogy(LoudnessTime, LoudnessShort, LoudnessTime, LoudnessDecayEnvelope, ':k', ...
        [LoudnessTime(Point1_TNU) LoudnessTime(Point2_TNU)], ...
        exp([BestFitLineTNU(Point1_TNU) BestFitLineTNU(Point2_TNU)]), '--xr');
    title(['TN(user): ', num2str(output.TN_User), ' s'])
    xlabel('Time (s)', 'fontsize', 15);
    ylabel('Loudness (sone)', 'fontsize', 15);
    %axis([0 3 1 4]);
    ymax = 10.^ceil(log10(max([max(LoudnessShort) max(BestFitLineTNU)])));
    ylim([ymax/1000 ymax])
end

warning ('on', 'MATLAB:Axes:NegativeDataInLogAxis');

if isstruct(RIR)
    doresultleaf(LoudnessShort,'Loudness [sone]',{'Time'},...
                 'Time',     LoudnessTime,       's',           true,...
                 'Function', {'Unique'}, 'categorical', [],...
                 'name','Loudness');
end

