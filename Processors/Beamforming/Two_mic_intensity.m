function [OUT varargout] = Two_mic_intensity(IN, spacing, temperature, atmosphericpressure, relativehumidity, calibrate, fftorder, flo, fhi,fade)
% Note: if cal field exist - currently cal field is returned for output,
% cal field for channel 2 (particle velocity) is meaningless and should be
% removed perhaps...
% The default mic calibration signal accounts for level
% differences, so only one gain cal field value is necessary. In this case
% cal field for chan 1 is appropriate for analysis.

% This function calibrates (if option is selected) a pair of pressure
% microphone signals and returns average sound pressure of the pair and
% particle velcocity vectors.
% Intensity and associated values (direction etc.) may be calculated with
% the sound intensity analyser.
% Input must be 2 channels.
% Default cal file: B&K type 4190 microphone with B&K type 2669
% preamplifier, fs = 48 kHz. freq range - 100-5 kHz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  MIC 1 %%%%%%%%%%%%%%%  MIC 2 %%%%%%%%
%%% B&K 4190 %% 2522008 %%%%%%%%%%%%%% 1811734 %%%%%%%
%%% B&K 2669 %% 2344491 %%%%%%%%%%%%%% 2417552 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% JH to double check this. Mics/pres possiblly still assembled.
% possibility of adding calibration option to record/calibrate and
% automaticically create IR's and applying appropiate fade in and fade out
% and add to audio3 field. streamline the process.
% test standing waves, test 2 speakers along orthogonal access in front and
% behind speakers. 
% 19mm spacing!!!!!!!!
%%
tic
if nargin ==1
    
    param = inputdlg({'Microphone Spacing (mm)';...
        'Temperature (deg C)';...
        'Atmospheric pressure (kPa)';...
        'Relative humidity (%)';...
        'Use calibration signal (0 (no cal) | 1 (user select) | 2 (default cal) ) ';... % option 2 load JH default cal
        'FFT order (power of 2)';... % not used
        'Low cut-off frequency';...
        'High cut-off frequency';...
        'N samples fade in/out'},... 
        'Settings',... % dialog window title.
        [1 50],...
        {'18';'20';'101.325';'50';'2';'16';'0';'20000';'512'}); % default values
    
    param = str2num(char(param));
    
    if length(param) < 9, param = []; end
    if ~isempty(param)
        spacing = param(1)/1000; % spacing in m
        temperature = param(2); % temperature in deg C
        atmosphericpressure = param(3)*1000; % pressure in pascals
        relativehumidity = param(4);
        calibrate = param(5);
        fftorder = param(6); %not implimented
        nfft = 2.^fftorder;
        flo = param(7);
        fhi = param(8);
        fade = param(9);
    end
else
    param = [];
end
%%
if isstruct(IN)
    fs = IN.fs;       % Extract the sampling frequency of the audio data
    IN = bandpass(IN, flo, fhi, 24, fs, 0); % Call AARAE's inbuilt function
                                            % bandpass filter. Default
                                            % settings are 24th order and
                                            % zero phase. check this.
    
    audio = IN.audio; % Extract the audio data
    
    
    %     % just an idea - not implemented at present
    %     if isfield(IN, 'audio3') && calibrate == 1
    %         cal_audio = IN.audio3;
    %     end
    %   %%
    
    if calibrate == 1 && ~exist('cal_audio','var')
        % USE CALIBRATION SIGNAL
        % 2-mic calibration signal (should be recorded with mic-pair cal
        % adaptor)
        
        % Use a menu & dialog box to select a wav file or audio within AARAE
        calfile = choose_audio; % call AARAE's choose_audio function
        if ~isempty(calfile)
            cal_audio = calfile.audio; % additional audio data
            fs2 = calfile.fs; % sampling rate
            
            if ~(fs2 == fs)
                disp('Unable to calibrate: calibration signal sampling rate does not match that of the input audio!')
            end
        end
    elseif calibrate == 2 && ~exist('cal_audio','var')
        calfile = importdata([cd '/Processors/Beamforming/Intensity1axis/intensity_cal_48k.mat']); % Default cal file.
        cal_audio = calfile.audio;
        fs2 = calfile.fs;
        if ~(fs2 == fs)
            disp('Unable to calibrate: calibration signal sampling rate does not match that of the input audio!')
        end
    end
    
    
elseif ~isempty(param) || nargin > 1
    audio = IN;
    fs = input_1;
end


if ~isempty(audio) && ~isempty(fs)
    
    [len,chans,bands] = size(audio); %take length here
    
    % ALLOW MULTIBAND?
    if bands > 1
        audio = sum(audio,3);
        disp('Multiband audio has been mixed into a single band')
    end
    %% apply half hann to start and end of input
    if  fade >0;
        fades = [hann(2*fade) hann(2*fade)];
        audio(1:fade,:) = audio(1:fade,:).* fades(1:fade,:); %fade in
        audio(len-fade+1:len,:) = audio(len-fade+1:len,:).*fades(fade+1:length(fades),:); % fade out
    else
    end
    
    % PHYSICAL PARAMETERS
    T = temperature + 273.15; % temperature in Kelvin
    psat = 6.1078 * 10.^(7.5*temperature/(temperature+237.3)); %saturation vapor pressure of water
    pv = relativehumidity * psat; % vapor pressure
    density = 1.2929 * (273.15 / T) * (atmosphericpressure - 0.3783*pv)/1.013e5;
    %         Weast, R.C. (Ed.) (1986) Handbook of Chemistry and Physics
    %         CRC Press Inc.; 67th edition; Boca Raton, Florida, USA.
    disp(['Density of air: ',num2str(density),' Kg/m^3'])
    
    
    % CALIBRATION
    if exist('cal_audio','var') && calibrate 
        [len2, chans2, bands2] = size(cal_audio); % new wave dimensions
        if chans2 == 2 && bands2 == 1
            totlen = len + len2 -1;
              fft1 = fft(cal_audio(:,1), totlen);
              fft2 = fft(cal_audio(:,2), totlen);
%             padlen = len - len2;
%             padlen = zeros(padlen,2);
%             cal_audio = [cal_audio; padlen];
%             fft1 = fft(cal_audio(:,1));
%             fft2 = fft(cal_audio(:,2));
            

            % Transfer function: crude, fast method. Actually more robust
            % than other methods which zero spectrum components and bugger
            % the phase response.
            % tf_cal12 = fft2./fft1;  % redundant
            tf_cal11 = fft1./fft1; % redundant, useful for conv().
            tf_cal21 = fft1./fft2;
            
            % apply calibration
            audio = [ifft(fft(audio(:,1,:),totlen).*tf_cal11) ...
                       ifft(fft(audio(:,2,:),totlen).*tf_cal21)]; %prep for time domain processing
            

        else
            disp('Unable to calibrate: calibration signal MUST be 2 channels, 1 band')
            
        end
        
    end
    
    
    % calculate acoustic quantities
    pressure = mean(audio,2);
    pressgradiant = (audio(:,2,:) - audio(:,1,:)) ./ spacing;
    velocity = -1/density * cumsum(pressgradiant)/fs;
    %%%%%%% Intensity vectors (and direction etc) to be calculated in analyser

    audio = [pressure, velocity];
    toc
    if isstruct(IN)
        OUT = IN; % You can replicate the input structure for your output
        OUT.audio = audio;
        OUT.properties.units = 'Pa | m/s';
        OUT.properties.units_type = 1;
        % The following line does not work in AARAE Release 8 because only 
        % one units_ref value is supported. It does work in AARAE Release 7
        % (because units_ref is not used in the GUI displays), and it works
        % in AARAE Release 9 (forthcoming).
        OUT.properties.units_ref = [2e-5 5e-8]; 
        OUT.chanID = {'pressure'; 'velocity'};
        %         OUT.cal = IN.cal(:,1);
        OUT.funcallback.name = 'Two_mic_intensity.m';
        OUT.funcallback.inarg = {spacing, temperature, atmosphericpressure, relativehumidity, calibrate, fftorder, flo, fhi,fade};
    else
        % You may increase the functionality of your code by allowing the
        % output to be used as standalone and returning individual
        % arguments instead of a structure.
        OUT = audio;
    end
    varargout{1} = fs;
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
% Copyright (c) 2016, Jonothan Holmes, Frederico Lopes Pereira, Densil
%                     Cabrera
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