function [OUT,varargout] = SII_IR(RIR, fs, P, BG, T, Method, presentation, v_effort, BandImportance, doplot)

% This function calculates the speech intelligibility index (SII)
% as per ANSI S3.5 (1997).
%
% The modulation transfer function of an impulse response is used in the
% calculation, in a similar way to the speech transmission index.
%
% While the various band approaches in the standard are all implemented,
% only the octave band approach is currently available when the function
% is run via dialog box.
%
%
% Code by Doheon Lee and Densil Cabrera
% version 1.01 (15 November 2013)

% INPUT: Mandatory input variables are 'RIR', 'P' and 'BG'.

%   RIR: a room impulse response. If a RIR has two channels, only left
%       channel is chosen for the SII calculations.

%   P: Combined speech and noise level as per Clause 3.17;
%       in critical bans from 150 Hz to 8.5 kHz (21 bands, see Table 1); or
%       in equally-contributing critical-bands from 350 Hz to 5.8 kHz (17 bands, see Table 2);
%       in 1/3 octave bands from 160 Hz to 8 kHz (18 bands, see Table 3); or
%       in octave bands from 250 Hz to 8 kHz (6 bands, see Table 4); or

%   BG: Background noise level in dB: (Default is 0 dB, under assumption that there is no external noise).
%       A string: a size of string must be matched with a size of 'P'.
%       0 or a string of zeros: when there is no background noise.

%   METHOD: a method for determining the Equivalent Speech Spectrum Level
%      and the Equivalent Noise Spectrum Level (Default is 1).
%       1: as per Clause 5.2
%       2: as per Clause 5.3

%   PRESENTATION: default is 1
%       1: Monaural listening
%       2: Bianrual listening

%   T: Equivalent Hearing Threshold Level (default is 0 dB, under assumption
%     that subjects are normal-hearing listeners)

%   V_EFFORT (default is 1): Determine the vocal effort
%       1: Normal (Tables 1-4)
%       2: Raised (Tables 1-4)
%       3: Loud (Tables 1-4)
%       4: Shout (Tables 1-4)

%   BANDIMPORTANCE:Band Importance (default is 1, note that this option is
%           not available for the equally-contributing ctiritical-bands procedure)
%       1: Average speech level as specified (Tables 1-4)
%       2: Nonsense Syllable test (Annex B)
%       3: Phonetically-balanced words of the CID-W22 (Annex B)
%       4: Monosyllables of the NU6 test (Annex B)
%       5: Diagnostic Rhyme Test (DRT) material (Annex B)
%       6: Short passages of easy reading material (Annex B)
%       7: Monosyllables of the speech in the presence of noise (SPIN) test (Annex B)

%   DOPLOT: Default is 0 when the first input is not a structure; Default
%   is 1 when the first input is a structure
%       0: Disable plots
%       1: plot the Equivalent Speech Spectrum Levels and Equivalent Noise
%               Spectrum Levls

%   FS: Audio sampling rate in Hz (only required if the first argument does
%   not contain both the audio data and sampling rate).

% OUTPUT:
%   S: Speech Intelligibility Index Value
%   Sav: Speech Intelligiblity Index considering Visual Cues (Annex B.1)
%   E: Equivalent Speech Spectrum Level (Clause 3.11)
%   N: Equivalent NOise Spectrum Level (Clause 3.13)


% EXAMPLES OF THE OCTAVE BAND PROCEDURE.
% P=[50 60 30 70 60 45]; ;BG=[30 20 25 30 30 30]; T=[10 10 10 10 10 10]
% Example 1: SII_IR('RIR1.wav', P, BG)
% Example 2: SII_IR('RIR1.wav', P, BG, 1, 2, T, 1, 1, 1)


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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    prompt = {'Method Clause 5.2 or 5.3 (1|2):', ...
        'Presentation monaural/binaural (1|2):', ...
        'Vocal Effort (1|2|3|4):', ...
        'Band importance (1:7):',...
        '250 Hz band hearing threshold (dB):', ...
        '500 Hz band hearing threshold (dB):', ...
        '1 kHz band hearing threshold (dB):', ...
        '2 kHz band hearing threshold (dB):',...
        '4 kHz band hearing threshold (dB):', ...
        '8 kHz band hearing threshold (dB):'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'1','2','1','1','0','0','0','0','0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        OUT = [];
        return
    else
        Method = str2num(answer{1,1});
        presentation = str2num(answer{2,1});
        v_effort = str2num(answer{3,1});
        BandImportance = str2num(answer{4,1});
        T = zeros(1,6);
        for k = 1:6
            T(k) = str2double(answer{k+4,1});
        end
    end
    
    prompt = {'250 Hz band level (dB):', ...
        '500 Hz band level (dB):', ...
        '1 kHz band level (dB):', ...
        '2 kHz band level (dB):',...
        '4 kHz band level (dB):', ...
        '8 kHz band level (dB):'};
    dlg_title = 'Combined Level of Speech & Noise';
    num_lines = 1;
    def = {'62','59','53','47','41','35'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        OUT = [];
        return
    else
        P = zeros(1,6);
        for k = 1:6
            P(k) = str2double(answer{k,1});
        end
    end
    
    prompt = {'250 Hz band level (dB):', ...
        '500 Hz band level (dB):', ...
        '1 kHz band level (dB):', ...
        '2 kHz band level (dB):',...
        '4 kHz band level (dB):', ...
        '8 kHz band level (dB):'};
    dlg_title = 'Background Noise Level';
    num_lines = 1;
    def = {'25','20','15','10','5','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        OUT = [];
        return
    else
        BG = zeros(1,6);
        for k = 1:6
            BG(k) = str2double(answer{k,1});
        end
    end

    doplot = 1;
end

% SETTINGS FOR AARAE
if isstruct(RIR)
    RIR = choose_from_higher_dimensions(RIR,1,1); 
    data = RIR.audio;
    fs = RIR.fs;
elseif ischar(RIR)
    [data,fs]=audioread(char(RIR));
else
    % direct input of an impulse response
    data = RIR;
    if exist('fs', 'var') == 0;
        % get the audio sampling rate if it has not been provided
        prompt = {'Audio sampling rate (Hz):'};
        dlg_title = 'Sampling rate';
        num_lines = 1;
        def = {'48000'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if isempty(answer)
            OUT = [];
            return
        else        
            fs = str2double(answer{1,1});           
        end
    end
end


if ~isempty(data) && ~isempty(fs) && ~isempty(P) && ~isempty(BG) && ~isempty(T) && ~isempty(Method) && ~isempty(presentation) && ~isempty(v_effort) && ~isempty(BandImportance) && ~isempty(doplot)
    % SETTINGS FOR SII CALCULATIONS
    if length(size(data))>2, data = data(:,1); end
    [len, chan]=size(data);
    if chan > 1
        data=data(:,1); % Select the first channel
    end;
    Nyquist=fs/2;
    time=((1:len)-1)'./fs;
    mf=[0.5, 1, 1.5, 2, 3, 4, 6, 8, 16];

    if nargin < 3 && ~isstruct(RIR)
        error('RIR, combined speech & noise levels, and noise levels are mendatory inputs')
        return
    end

    if max(BG-P)>0
        error('Noise level cannot be greater than combined speech and noise level!')
        return
    end

    if exist('Method', 'var') == 0;
        Method=1; % Clause 5.2 by default
    end

    if exist('T', 'var') == 0;
        T=zeros(1,length(P)); %Threshold Level in dB HL for normal-hearing listeners.
    end

    if exist('presentation', 'var') == 0;
        presentation = 1;
    end

    if presentation == 2
        T=T-1.7; % Adjustment for the binaural listening (Clause 5.2.4)
    end

    if exist('v_effort', 'var') == 0;
        v_effort=1;
    end

    if exist('BandImportance', 'var') ==0;
        BandImportance=1;
    elseif BandImportance ~=1 && length(P) == 17
        warning('Band Importance is changed to 1, \n as the others are not available for the equally-conbributing critical-band procedure', 7);
        BandImportance=1;
    end

    if exist('doplot', 'var') == 0
        doplot=0;
    end


    switch length(P)
        case 21 % CRITICAL BAND PROCEDURE (Table 1)
            fc=[150, 250, 350, 450, 570, 700, 840, 1000, 1170, 1370, 1600, ...
                1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800, 7000, 8500];
            f_low=[100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, ...
                1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700];
            f_hi=[200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, ...
                2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500];
            H=[0.6, 1.0, 1.4, 1.4, 1.9, 2.8, 3.0, 2.6, 2.6, 3.6, 6.1, 10.5, ...
                13.8, 16.8, 15.8, 14.9, 14.3, 12.4, 7.9, 4.3, 0.5];
            X=[1.5, -3.9, -7.2, -8.9, -10.3, -11.4, -12.0, -12.5, -13.2, ...
                -14.0, -15.4, -16.9, -18.8, -21.2, -23.2, -24.9, -25.9, ...
                -24.2, -19.0, -11.7, -6.0];
            U=[31.44, 34.06, 34.21, 28.69; ...
                34.75, 38.98, 41.55, 42.50; ...
                34.14, 38.62, 43.68, 47.14; ...
                34.58, 39.84, 44.08, 48.46; ...
                33.17, 39.44, 45.34, 50.17; ...
                30.64, 37.99, 45.22, 51.68; ...
                27.59, 35.85, 43.60, 51.43; ...
                25.01, 33.86, 42.16, 51.31; ...
                23.52, 32.56, 41.07, 49.40; ...
                22.28, 30.91, 39.68, 49.03; ...
                20.15, 28.58, 37.70, 47.65; ...
                18.29, 26.37, 35.62, 45.47; ...
                16.37, 24.34, 33.17, 43.13; ...
                13.80, 22.35, 30.98, 40.80; ...
                12.21, 21.04, 29.01, 39.15; ...
                11.09, 19.56, 27.71, 37.30; ...
                9.33, 16.78, 25.41, 34.41; ...
                5.84, 12.14, 19.20, 29.01; ...
                3.47,  9.04, 15.37, 25.17; ...
                1.78,  6.36, 12.61, 22.08; ...
                -0.14,  3.44,  9.62, 18.76];
            I=[0.0103, 0.0000, 0.0507, 0.0234, 0.0122, 0.0192, 0.0130; ...
                0.0261, 0.0230, 0.0677, 0.0368, 0.0553, 0.0312, 0.0478; ...
                0.0419, 0.0385, 0.0641, 0.0520, 0.0581, 0.0926, 0.0451; ...
                0.0577, 0.0410, 0.0552, 0.0672, 0.0672, 0.1031, 0.0470; ...
                0.0577, 0.0433, 0.0474, 0.0638, 0.0680, 0.0735, 0.0523; ...
                0.0577, 0.0472, 0.0468, 0.0566, 0.0667, 0.0611, 0.0591; ...
                0.0577, 0.0473, 0.0466, 0.0503, 0.0587, 0.0495, 0.0591; ...
                0.0577, 0.0470, 0.0502, 0.0465, 0.0547, 0.0440, 0.0503; ...
                0.0577, 0.0517, 0.0586, 0.0539, 0.0563, 0.0440, 0.0503; ...
                0.0577, 0.0537, 0.0591, 0.0576, 0.0575, 0.0490, 0.0556; ...
                0.0577, 0.0582, 0.0586, 0.0642, 0.0625, 0.0486, 0.0699; ...
                0.0577, 0.0679, 0.0609, 0.0741, 0.0598, 0.0493, 0.0625; ...
                0.0577, 0.0745, 0.0596, 0.0709, 0.0555, 0.0490, 0.0602; ...
                0.0577, 0.0750, 0.0618, 0.0621, 0.0521, 0.0547, 0.0684; ...
                0.0577, 0.0685, 0.0501, 0.0553, 0.0480, 0.0555, 0.0638; ...
                0.0577, 0.0662, 0.0439, 0.0505, 0.0443, 0.0493, 0.0605; ...
                0.0577, 0.0636, 0.0370, 0.0417, 0.0356, 0.0359, 0.0534; ...
                0.0460, 0.0607, 0.0268, 0.0291, 0.0280, 0.0387, 0.0394; ...
                0.0343, 0.0511, 0.0201, 0.0186, 0.0237, 0.0256, 0.0291; ...
                0.0226, 0.0216, 0.0186, 0.0141, 0.0205, 0.0219, 0.0132; ...
                0.0110, 0.0000, 0.0612, 0.0113, 0.0153, 0.0043, 0.0000];

        case 17 %EQUALLY-CONTRIBUTIG CRIBICAL-BAND PROCEDURE (TABLE 2)
            fc=[350, 450, 570, 700, 840, 1000, 1170, 1370, 1600, 1850, ...
                2150, 2500, 2900, 3400, 4000, 4800, 5800];
            f_low=[300, 400, 510, 630, 770, 920, 1080, 1270, 1480, ...
                1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300];
            f_hi=[400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, ...
                2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400];
            H=[1.4, 1.4, 1.9, 2.8, 3.0, 2.6, 2.6, 3.6, 6.1, 10.5, 13.8, ...
                16.8, 15.8, 14.9, 14.3, 12.4, 7.9];
            X=[-7.2, -8.9, -10.3, -11.4, -12.0, -12.5, -13.2, ...
                -14.0, -15.4, -16.9, -18.8, -21.2, -23.2, -24.9, -25.9, ...
                -24.2, -19.0];
            U=[34.14, 38.62, 43.68, 47.14; ...
                34.58, 39.84, 44.08, 48.46; ...
                33.17, 39.44, 45.34, 50.17; ...
                30.64, 37.99, 45.22, 51.68; ...
                27.59, 35.85, 43.60, 51.43; ...
                25.01, 33.86, 42.16, 51.31; ...
                23.52, 32.56, 41.07, 49.40; ...
                22.28, 30.91, 39.68, 49.03; ...
                20.15, 28.58, 37.70, 47.65; ...
                18.29, 26.37, 35.62, 45.47; ...
                16.37, 24.34, 33.17, 43.13; ...
                13.80, 22.35, 30.98, 40.80; ...
                12.21, 21.04, 29.01, 39.15; ...
                11.09, 19.56, 27.71, 37.30; ...
                9.33, 16.78, 25.41, 34.41; ...
                5.84, 12.14, 19.20, 29.01; ...
                3.47,  9.04, 15.37, 25.17];
            I=repmat(0.0588, 17,1);

        case 18 % 1/3 OCTAVE BAND PROCEDURE (TABLE 3)
            bandnumber=22:1:39;
            fc=10.^(bandnumber/10);
            f_low=fc./10^(0.15/3);
            f_hi=fc.*10^(0.15/3);
            H=[0, 0.5, 1.0, 1.4, 1.5, 1.8, 2.4, 3.1, 2.6, 3.0, 6.1, 12.0, 16.8, ...
                15.0, 14.3, 10.7, 6.4, 1.8];
            X=[0.6, -1.7, -3.9, -6.1, -8.2, -9.7, -10.8, -11.9, -12.5, -13.5, ...
                -15.4, -17.7, -21.2, -24.2, -25.9, -23.6, -15.8, -7.1];
            U=[32.41	33.81	35.29	30.77;
                34.48	33.92	37.76	36.65;
                34.75	38.98	41.55	42.50;
                33.98	38.57	43.78	46.51;
                34.59	39.11	43.30	47.40;
                34.27	40.15	44.85	49.24;
                32.06	38.78	45.55	51.21;
                28.30	36.37	44.05	51.44;
                25.01	33.86	42.16	51.31;
                23.00	31.89	40.53	49.63;
                20.15	28.58	37.70	47.65;
                17.32	25.32	34.39	44.32;
                13.18	22.35	30.98	40.80;
                11.55	20.15	28.21	38.13;
                9.33	16.78	25.41	34.41;
                5.31	11.47	18.35	28.24;
                2.59	7.67	13.87	23.45;
                1.13	5.07	11.39	20.72];
            I=[0.0083	0.0000  0.0365	0.0168	0.0000	0.0114	0.0000; ...
                0.0095	0.0000	0.0279	0.0130	0.0240	0.0153	0.0255; ...
                0.0150	0.0153	0.0405	0.0211	0.0330	0.0179	0.0256; ...
                0.0289	0.0284	0.0500	0.0344	0.0390	0.0558	0.0360; ...
                0.0440	0.0363	0.0530	0.0517	0.0571	0.0898	0.0362; ...
                0.0578	0.0422	0.0518	0.0737	0.0691	0.0944	0.0514; ...
                0.0653	0.0509	0.0514	0.0658	0.0781	0.0709	0.0616; ...
                0.0711	0.0584	0.0575	0.0644	0.0751	0.0660	0.0770; ...
                0.0818	0.0667	0.0717	0.0664	0.0781	0.0628	0.0718; ...
                0.0844	0.0774	0.0873	0.0802	0.0811	0.0672	0.0718; ...
                0.0882	0.0893	0.0902	0.0987	0.0961	0.0747	0.1075; ...
                0.0898	0.1104	0.0938	0.1171	0.0901	0.0755	0.0921; ...
                0.0868	0.1120	0.0928	0.0932	0.0781	0.0820	0.1026; ...
                0.0844	0.0981	0.0678	0.0783	0.0691	0.0808	0.0922; ...
                0.0771	0.0867	0.0498	0.0562	0.0480	0.0483	0.0719; ...
                0.0527	0.0728	0.0312	0.0337	0.0330	0.0453	0.0461; ...
                0.0364	0.0551	0.0215	0.0177	0.0270	0.0274	0.0306; ...
                0.0185	0.0000	0.0253	0.0176	0.0240	0.0145	0.0000];

        case 6 %OCTAVE-BAND PROCEDURE (TABLE 4)
            bandnumber=24:3:39;
            fc=10.^(bandnumber/10);
            f_low=fc./10^0.15;
            f_hi=fc.*10^0.15;
            H=[1.0, 1.8, 2.6, 12.0, 14.3, 1.8];
            X=[-3.9, -9.7, -12.5, -17.7, -25.9, -7.1];
            U=[34.75 38.98 41.55 42.50;
                34.27 40.15 44.85 49.24;
                25.01 33.86 42.16 51.31;
                17.32 25.32 34.39 44.32;
                9.33  16.78 25.41 34.41;
                1.13   5.07 11.39 20.72];
            I=[0.0617,0.0437, 0.1549, 0.0853, 0.0960, 0.1004, 0.0871; ...
                0.1671, 0.1294, 0.1562, 0.1912, 0.2043, 0.2551, 0.1493; ...
                0.2373, 0.2025, 0.2165, 0.2110, 0.2343, 0.1960, 0.2206; ...
                0.2648, 0.3117, 0.2768, 0.3090, 0.2643, 0.2322, 0.3022; ...
                0.2142, 0.2576, 0.1488, 0.1682, 0.1501, 0.1744, 0.2102; ...
                0.0549, 0.0551, 0.0468, 0.0353, 0.0510, 0.0419, 0.0306];
    end

    % ONE-THIRD, ONE OCTAVE OR CRITICAL BAND FILTER
    P_octave=zeros(len, length(fc));
    if length(fc) == 6
        for i=1:length(fc);
            warning('off');
            H=fdesign.octave(1, 'Class 0', 'N,F0', 6, fc(i), fs);
            Hd=design(H, 'butter', 'SOSScaleNorm', 'Linf');
            warning('on');
            P_octave(:,i)=filter(Hd, data);
        end
    elseif length(fc) == 18
        for i=1:length(fc);
            warning('off');
            H=fdesign.octave(3, 'Class 0', 'N,F0', 8, fc(i), fs);
            Hd=design(H, 'butter', 'SOSScaleNorm', 'Linf');
            warning('on');
            P_octave(:,i)=filter(Hd, data);
        end
    else
        for i=1:length(fc);
            wn=[f_low(i)/Nyquist, f_hi(i)/Nyquist];
            [b, a]=butter(3, wn);
            P_octave(:,i)=filtfilt(b, a, data);
        end
    end

    % MODULATION TRANSFER FUNCTION
    MTF=zeros(length(mf), length(fc));
    for k=1:length(fc);
        for j=1:length(mf)
            MTF_num=abs(sum(P_octave(:,k).^2.*exp(-2i*pi*mf(j).*time)));
            MTF_den=sum(P_octave(:,k).^2);
            MTF(j,k)=MTF_num/MTF_den;
        end
    end


    % APPARENET SPEECH TO NOISE RATIO
    R=10*log10((MTF+eps)./(1-MTF+eps));
    R(R<-15)=-15; R(R>15)=15;
    R=mean(R);

    % EQUAVALENT SPEECH AND NOISE SPECTRUM
    switch Method
        case 1
            E=R+10*log10( 10.^(0.1*P)./(1+10.^(0.1*R)));
            N=10*log10(10.^((E-R)/10)+10.^(BG/10)); % Note in Clause 5.2.3
        case 2
            E=R+10*log10( 10.^(0.1*P)./(1+10.^(0.1*R)))-H;
            N=10*log10(10.^((E-R)/10)+10.^(BG/10))-H;
    end

    % if exist('BG', 'var') == 0 || isequal(BG, zeros(1, length(P)));
    %     N=repmat(-50, 1, length(P)); % Note in Clause 4.2
    % end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% SII Calculation %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Z=N; % EQUIVALENT MASKING SPECTRUM
    V=E-24; % SELF-SPEECH MASKING SPECTRUM
    B=max(Z, V); % SELECT LARGER OF Z OR V;
    Z(1)=B(1);

    Z1=zeros(1,length(fc)-1);
    Z2=zeros(1, length(fc)-1);
    switch length(P)
        case {17 21}
            C=-80+0.6*(B+10*log10(f_hi-f_low)); % THE SLOPE PER OCTAVE OF THE SPREAD OF MASKING
            for i=2:length(fc);
                Z1(i)=10^(0.1*N(i));
                Z2(i)=sum(10.^(0.1*(B(1:i-1) +3.32.*C(1:i-1) .* log10(fc(i)./f_hi(1:i-1)))));
                Z(i)=10*log10(Z1(i)+Z2(i)); % EQUIVALENT MASKING SPECTRUM LEVEL
            end
        case {18}
            C=-80+0.6*(B+10*log10(fc)-6.353); % THE SLOPE PER OCTAVE OF THE SPREAD OF MASKING
            for i=2:length(fc);
                Z1(i)=10^(0.1*N(i));
                Z2(i)=sum(10.^(0.1*(B(1:i-1) +3.32.*C(1:i-1) .* log10(0.89*fc(i)./fc(1:i-1)))));
                Z(i)=10*log10(Z1(i)+Z2(i)); % EQUIVALENT MASKING SPECTRUM LEVEL
            end
        case {6}
            % FOR THE OCTAVE BAND PROCEDURE, THE EQUIVALENET MASKING SPECTRUM MUST BE
            % EQUAL TO THE EQUIVALENT NOISE SPECTRUM (ANNEX C).
            if length(P)==6
                Z=N; % EQUIVALENT MASKING SPECTRUM LEVEL
            end
    end


    X=X+T;% EQUIVALENT INTERNAL NOISE SPECTRUM LEVEL
    D=max(Z, X); % EQUIVALENT DISTURBANCE SPECTRUM LEVEL
    L=1-(E-U(:,v_effort)'-10)/160; L(L>1)=1; %LEVEL DISTORTION FACTOR
    K=(E-D+15)/30; K(K>1)=1; K(K<0)=0; %TEMPORARY VARIABLE
    A=L.*K; %BAND AUDIBILITY FUNCTION
    S=sum(I(:,BandImportance)'.*A); %SPEECH INTELLIGIBLITY INDEX

    % SII WITH VISUAL CUES (ANNEX B.1)
    if S <= 0.2
        Sav=0.1+1.5*S;
    elseif S > 0.2
        Sav=0.25+0.75*S;
    end
    
    if isstruct(RIR)
        OUT.S = S;
        OUT.Sav = Sav;
        OUT.E = E;
        OUT.N = N;
        OUT.funcallback.name = 'SII_IR.m';
        OUT.funcallback.inarg = {fs,P,BG,T,Method,presentation,v_effort,BandImportance,doplot};
    else
        OUT  = [];
        varargout{1} = S;
        varargout{2} = Sav;
        varargout{3} = E;
        varargout{4} = N;
    end
    
    if isstruct(RIR)
        doresultleaf(MTF,'Modulation TF',{'Modulation_frequency'},...
                     'Modulation_frequency', num2cell(mf),                            'Hz', true,...
                     'Frequency',            num2cell([250,500,1000,2000,4000,8000]), 'Hz', false,...
                     'name','Modulation_TF');

        doresultleaf([P',E',N',BG'],'SPL [dB]',{'Frequency'},...
                     'Frequency', num2cell(fc),                                              'Hz',          true,...
                     'Level',     {'Received SPL','Eq Spch Lvl','Eq Noise Lvl','Noise Lvl'}, 'categorical', [],...
                     'name','Band_SPL');
    end
                 
    % PLOTTING
    if doplot
        figure('Name', ['Speech Intelligibility Index: ', num2str(S)])


        % plot of mtfs
        if length(P) == 6

            subplot(2,1,1)

            % rainbow colours
            colors = [255, 0, 0; ... % red
                255, 128, 0; ... % orange
                204, 204, 0; ... % dark yellow
                0, 204, 0; ... % mid green
                0, 204, 204; ... % dark cyan
                0, 0, 255; ... % blue
                127, 0, 255]; % violet
            colors = colors / 255; % rescale to 0-1 range

            hold on
            plot(mf,MTF(:,6),'Color',colors(7,:),'DisplayName','8 kHz');
            plot(mf,MTF(:,5),'Color',colors(6,:),'DisplayName','4 kHz');
            plot(mf,MTF(:,4),'Color',colors(5,:),'DisplayName','2 kHz');
            plot(mf,MTF(:,3),'Color',colors(4,:),'DisplayName','1 kHz');
            plot(mf,MTF(:,2),'Color',colors(3,:),'DisplayName','500 Hz');
            plot(mf,MTF(:,1),'Color',colors(2,:),'DisplayName','250 Hz');

            set(gca, 'Xscale','log')
            xlim([0.5 16])
            ylim([0 1])

            % legend
            legend('show','Location','EastOutside');

            % Create xlabel
            xlabel('Modulation Frequency (Hz)');

            title('Modulation Transfer Function');
            grid on
            hold off
        end

        if length(P) == 6
            subplot(2,1,2)
        end
        hold on

        % Sound pressure level of the received signal (& noise)
        plot(fc,P,'b','Marker','o','DisplayName','Received')

        % Equivalent speech level
        plot(fc,E,'r','Marker','+','DisplayName','Eq Spch Lvl')

        % Equivalent noise level
        plot(fc,N,'Color', [0 0.6 0],'Marker','x','DisplayName','Eq Noise Lvl')

        % Noise level
        plot(fc,BG,'k','Marker','x','DisplayName','Noise Lvl')

        set(gca, 'Xscale','log')
        xlim([100 10000])

        legend('show','Location','EastOutside');

        ylabel('SPL (dB)')
        xlabel('Frequency (Hz)');
        title('Band Sound Pressure Level')
        grid on
        hold off

        if length(P) == 6

            f = figure('Name',['Speech Intelligibility Index: ', num2str(S)], ...
                'Position',[200 200 620 360]);
            %[left bottom width height]
            dat1 = MTF;
            cnames1 = {'250', '500', '1k', '2k', '4k', '8k'};
            rnames1 = {'0.5', '1', '1.5', '2', '3', '4', '6', '8', '16'};
            t1 =uitable('Data',dat1,'ColumnName',cnames1,'RowName',rnames1);
            %set(t,'ColumnWidth',{60});

            dat2 = [P;E;N;BG];
            cnames2 = {'250', '500', '1k', '2k', '4k', '8k'};
            rnames2 = {'Received SPL (dB)','Eq Speech Level (dB)', ...
                'Eq Noise Level (dB)', 'Noise Level (dB)'};
            t2 =uitable('Data',dat2,'ColumnName',cnames2,'RowName',rnames2);
            %set(t,'ColumnWidth',{60});

            [~,tables] = disptables(f,[t2 t1],{'Sound levels','SII'});
            OUT.tables = tables;
        end
    end
end
