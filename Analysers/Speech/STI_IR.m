function [Verbose,varargout]=STI_IR(IR, fs, Lsignal, Lnoise, AuditoryMasking, NoiseCorrection, doplot, FilterVersion,FilterStrength,FilterComp)
% This function calculates the speech transmission index (STI)
% in accordance with IEC 60268-16 (2011-06).
% It implements the 'indirect method of measuring STI using the impulse
% response' (section 6 of IEC 60268-16 (2011)).
%
% The impulse response (IR) (of the room or system) must be made prior to
% calling this function, in accordance with ISO 18233. The IR should be at
% least 1.6 s long, and not shorter than half the reverberation time of the
% room (in each octave band). The IR should be noise-free in each octave
% band (or at least SNR >= 20 dB). It is advisable to check that the IR is
% correctly prepared prior to analysis.
%
% The signal and noise octave band sound pressure levels also need to be
% determined as described in IEC 60268-16 (2011),
% unless it is assumed that noise is negligible.
%
% Reference: D. Cabrera, D. Lee, G. Leembruggen & D. Jimenez (in press
% 2014), "Increasing robustness in the calculation of the speech transmission
% index from impulse responses", Building Acoustics 21(3).
%
% Code by Doheon Lee and Densil Cabrera
% version 1.10 (22 July 2014)
%
% * For further information type 'doc STI_IR' in the MATLAB command line and
%   click on 'View code for STI_IR'.

% INPUT ARGUMENTS
%
%   If no input arguments are used in the function call, then
%   an impulse response wave file can be opened via a dialog box, and further
%   function settings are also entered via dialog boxes.
%
% IR: an impulse response, which can be a wave file name, or
%   data (vector  or matrix) or a structure.
%   If it is data, then fs needs to be specified if it is not
%   48 kHz (input argument 7). Multichannel IRs are analysed
%   sequentially (in the order of the channels).
%   If it is a structure, then the structure should contain at least the
%   following leaves:
%     IR.audio, which is the impulse response wave data.
%     IR.fs, which is the audio sampling rate in Hz.
%   If the input is a structure, then all other input arguments are
%   ignored, and user settings are entered using a sequence of dialog
%   boxes. The purpose of this is so that the function can be used directly
%   within a larger project with a GUI.
%   This argument can alternatively be a special input code, 'v', which is
%   used to run the function in valdation mode. In this mode, data from Annex
%   M of the standard (a noise-free MTF, and the operational speech and noise
%   levels) are used to compare the processing with the results in Annex M.
%   Discrepancies appear to be due to precision of the data in the Annex M
%   table.
%
% Lsignal: the sound pressure level of speech in octave bands from
%   125 Hz to 8 kHz. Please check Annex J of IEC 60268-16 (2011) to
%   correctly measure and/or prepare this data.
%
% Lnoise: the sound pressure level of background noise in octave bands from
%   125 Hz to 8 kHz.
%   When a multichannel IR is input, Lsignal and/or Lnoise can either:
%    (i) have the same number of rows as input channels (where the 7
%    columns are used for the 7 octave bands from 125 Hz to 8 kHz); or
%    (ii) have a single row of values, in which case these are used to
%    derive STI for all of the IR channels.
%
% AuditoryMasking
%   0: disable the auditory masking functions
%   1 (or 2011): enable the auditory masking functions (default), following
%     IEC 60268-16 (2011).
%   2 (or 2003): enable the auditory masking functions, following
%     IEC 60268-16 (2003).
%   -2011: enable auditory masking from IEC 60268-16 (2011), but without
%     auditory reception threshold
%   -2003: enable auditory masking from IEC 60268-16 (2003), but without
%     auditory reception threshold.
%   -1: enable auditory reception threshold, but without auditory masking.
%   1985: no masking or auditory reception threshold, with S/N valuese
%     calculated as described by Houtgast (1985).
%
% NoiseCorrection:
%   0: disable adjustment for background noise
%   1: enable the adjustment for background noise (default) - i.e. adjustment
%   of the modulation transfer function for the specified background noise
%   levels in relation to the specified speech levels.
%   Disabling NoiseCorrection should yield almost the same result as for a
%   very low background noise level (via the 'Lnoise' input argument).
%
% doplot
%   0: do not plot; do output AARAE result leafs
%   1: plot the modulation transfer functions and transmission
%   indices(default), and do output AARAE result leafs
%   -1: do plot but do not output AARAE result leafs
%   -2: do not plot and do not output AARAE result leafs
%
% fs: the audio sampling rate of the impulse response.
%   If a wave file is read, then this is not used (fs is read from the file).
%   If an IR is input directly, then fs specifies its sampling rate, but if
%   omitted, then fs defaults to 48 kHz. The audio sampling rate is recommended
%   to be at least 44.1 kHz (filtering will fail if the sampling rate is
%   too low). If a structure is input in argument 1, then fs must be
%   specified within that structure.
%
% FilterVersion
%   0: use 8x 1st order Butterworth octave band filters, generated using
%      explicit simple code with filtfilt to make them linear phase 
%   1: use AARAE's 6th order Butterworth octave band filters, generated using
%      Matlab's filterbuilder (these filters are not linear phase).
%   2: use AARAE's linear phase filters (this is the method that best
%      complies with the standard). The filter pseudo-order can be edited 
%      in the code (default)
%   Note that the differences between these approaches can be seen as small
%   deviations in the MTF values, but they are unlikely to affect
%   calculated STI values.
%
% FilterStrength
%   Applies factor to the filter order (or pseudo-order) - i.e. a value of
%   1 yields the filter orders described above, 2 doubles the filter order.
%   Using the current settings in the code, a value of 1 yields the
%   equivalent of 12th order filterbank response (72 dB/octave), a value of
%   2 yields 24th order (144 dB/octave), etc.
%   THIS ONLY APPLIES TO FilterVersion = 2
%
% FilterComp
%   0: no filter compensation
%   1: calculates the MTF of the filterbank itself using a delta function
%   within silence of the same period as the impulse response being tested.
%   (The delta function is in the middle of the silent period).
%   Then the MTF of the impulse response is divided by the MTF of the
%   filterbank, so that the filterbank's own MTF is removed from the
%   result. This probably will not significantly affect the STI values, but
%   it may improve the accuracy of the MTF itself (if high precision is
%   required). Filter compensation is most useful if high order filters are
%   used.
%   THIS ONLY APPLIES TO FilterVersion = 2
%
% OUTPUT ARGUMENTS
% M_STI is the male speech transmission index. Values are in order of the
%   IR channels.
%
% F_STI is the female speech transmission index.  Values are in order of the
%   IR channels.
%
% Verbose is an output structure containing male and female STI and
%   STIPA(IR) values, octave-band modulation transfer indices,
%   and the modulation transfer functions from which they
%   were derived. The suffix '0' indicates values for which no background
%   noise correction nor auditory masking treatment was applied.
%   The seven columns of the modulation transfer functions
%   represent the octave band carrier frequencies from 125 Hz to 8 kHz,
%   i.e. [125, 250, 500, 1000, 2000, 4000, 8000].
%   The 14 rows represent the 1/3-octave spaced modulation frequencies
%   from 0.63 Hz to 12.5 Hz, i.e.
%   [0.63; 0.8; 1; 1.25; 1.6; 2; 2.5; 3.15; 4; 5; 6.3; 8; 10; 12.5].
%   The third dimension is used for multi-channel IRs,
%   in order of the channels.
%
% The function also outputs a set of charts (if doplot is set), and
% generates a simple auralization (if set by the dialog box). The
% auralization consists of a user-selected wav file convolved with the
% impulse response (1 or 2 channels), which is then filtered to match
% Lsignal and mixed with noise that has been filtered to match Lnoise.
%
%
% EXAMPLES OF CALLING THIS FUNCTION
%
%   To run the function only using dialog box data entry:
% [M_STI,F_STI,Verbose]=STI_IR;
%
%   If RIR_room1.wav is a wave file containing an impulse response (only),
%   and is in the current directory, then the following could be used
%   (with signal and noise levels as specified):
% [M_STI,F_STI,Verbose]=STI_IR('RIR_room1.wav', [60 70 70 80 60 30 20], [50 50 60 65 30 25 20]);
%
%   If the IR is already in the Matlab Workspace (and is called 'data'),
%   then the following could be used:
% [M_STI,F_STI]=STI_IR(data, [60 70 70 80 60 30 20], [50 50 60 65 30 25 20]);
%   (This assumes fs = 48 kHz, which is the default sampling rate.)
%
%   If the IR is already in Matlab, but has a sampling rate of 44.1 kHz, then
%   almost all of the input arguments must be used, e.g.:
% [M_STI,F_STI]=STI_IR(data, [60 70 70 80 60 30 20], [50 50 60 65 30 25 20],1,1,1,44100);
%   (In the above, auditory masking and noise correction are employed, and
%   plotting is done.)
%
%   For a multichannel IR, corresponding multichannel signal and noise
%   values can be specified, e.g., for a 2 channel IR:
%  S = [60 70 70 80 60 30 20; 50 60 60 70 50 20 10];
%  N = [70 65 60 55 50 45 40; 65 60 55 50 45 40 35];
%  [M_STI,F_STI,Verbose]=STI_IR(data,S,N);
%   If only one set of signal and noise levels is supplied, then it will be
%   applied to each channel of the multichannel IR.
%
%   To call the function with an IR and fs in a structure, but with other
%   settings controlled via dialog boxes:
%  [M_STI,F_STI,Verbose]=STI_IR(IR);
%   where IR.audio contains the impulse response waveform, and IR.fs
%   contains the audio sampling rate.
%
%   To run the function in validation mode:
%  [M_STI,F_STI,Verbose]=STI_IR('v');


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

%***************************************************************
% Check the input arguments and data

% The dialogmode flag indicates whether dialog boxes will be used for inputs
% It is cleared here to 0 (false) to simplify the subsequent code, but if
% there are no input arguments or if the first input argument
% is a structure, then it is set to 1 (true) when
% that is checked (later in the code).
dialogmode = 0;

% The validate flag indicates whether validation data from Annex M will be
% tested. It is set later if the first input argument is 'v'.
validate = 0;

% Flag for auralization, which can be set in dialogmode.
doAuralization = 0;

Verbose.tables = [];


if nargin == 0
    % DIALOG BOX TO IMPORT AUDIO DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Use a dialog box to select a wav file.
    % The variable IR is used to store the file name.
    [IR, pathname] = uigetfile({'*.wav'}, ...
        'Open a wav file containing an impulse response');
    
    % Open the wav file, and acquire IR data and sampling rate.
    [data, fs] = audioread([pathname,IR]);
    
    % If you have an older version of Matlab that does not support audioread,
    % then use the following instead:
    %[data, fs] = wavread([pathname,IR]);
    
    % Set dialogmode.
    dialogmode = 1;
end

if (isstruct(IR) || dialogmode) && nargin < 2
    % DIALOG BOX ENTRY MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If a structure is input or if dialogmode has been set, then
    % this function runs with dialog boxes to input the main parameters
    % (but not the IR waveform and fs, which need to be supplied
    % by the structure if they have not already been imported).
    % The structure should (at least) contain the following two leaves:
    %    IR.audio is the impulse response
    %    IR.fs is the audio sampling rate in Hz
    
    dialogmode = 1;
    doplot = 1;
    
    % required data from input structure
    if isstruct(IR)
        IR = choose_from_higher_dimensions(IR,2,1); 
        data = IR.audio;
        fs = IR.fs;
    end
    
    % dialog box to get auditory masking, noise correction  and
    % auralization settings
    prompt = {'Auditory Masking (0 | 1 | 2):', ...
        'Noise correction (0 | 1):', ...
        'Filter Strength Factor:', ...
        'Compensate for filterbank''s own MTF (0 | 1):', ...
        'Auralization (0 | 1):'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'1','1','1','0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(answer)
        M_STI = [];
        F_STI = [];
        Verbose = [];
        return
    else
        AuditoryMasking = str2num(answer{1,1});
        NoiseCorrection = str2num(answer{2,1});
        FilterStrength = str2num(answer{3,1});
        FilterComp = str2num(answer{4,1});
        doAuralization = str2num(answer{5,1});
    end
    
    
    % dialog box to get signal level in octave bands
    if ~AuditoryMasking == 0 || NoiseCorrection
        
        prompt = {'125 Hz band level (dB):', ...
            '250 Hz band level (dB):', ...
            '500 Hz band level (dB):', ...
            '1 kHz band level (dB):', ...
            '2 kHz band level (dB):',...
            '4 kHz band level (dB):', ...
            '8 kHz band level (dB):'};
        dlg_title = 'Speech';
        num_lines = 1;
        def = {'62.9','62.9','59.2','53.2','47.2','41.2','35.2'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if isempty(answer)
            M_STI = [];
            F_STI = [];
            Verbose = [];
            return
        else
            Lsignal = zeros(1,7);
            for k = 1:7
                Lsignal(k) = str2double(answer{k,1});
            end
        end
    else
        Lsignal = 60 + [2.9, 2.9, -0.8, -6.8, -12.8, -18.8, -24.8];
    end
    
    
    % dialog box to get noise level in octave bands
    if NoiseCorrection
        
        prompt = {'125 Hz band level (dB):', ...
            '250 Hz band level (dB):', ...
            '500 Hz band level (dB):', ...
            '1 kHz band level (dB):', ...
            '2 kHz band level (dB):', ...
            '4 kHz band level (dB):', ...
            '8 kHz band level (dB):'};
        dlg_title = 'Noise';
        num_lines = 1;
        def = {'50','45','40','35','30','25','20'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if isempty(answer)
            M_STI = [];
            F_STI = [];
            Verbose = [];
            return
        else
            Lnoise = zeros(1,7);
            for k = 1:7
                Lnoise(k) = str2double(answer{k,1});
            end
        end
    else
        Lnoise = [-10 -10 -10 -10 -10 -10 -10];
        
    end
    
elseif ischar(IR)
    if IR == 'v'
        % VALIDATION MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % in this mode, the noise-free modulation transfer function, and
        % operational signal and noise levels from Annex M are used. (No impulse
        % response is analysed.)
        validate = 1;
        
        % Data from Annex M
        
        Lsignal = [82.9 82.9 79.2 73.2 67.2 61.2 55.2];
        Lnoise = [55.5 47.5 41.5 37.5 34.5 32.5 30.5];
        
        % noise-free mtf
        MTF_ch = [0.983 0.960 0.978 0.990 0.990 0.986 0.997; ... % 0.63 Hz
            0.968 0.936 0.959 0.974 0.980 0.979 0.995; ... % 0.8 Hz
            0.947 0.904 0.931 0.953 0.966 0.968 0.992; ... % 1 Hz
            0.920 0.869 0.898 0.927 0.949 0.955 0.987; ... % 1.25 Hz
            0.886 0.826 0.852 0.892 0.925 0.935 0.981; ... % 1.6 Hz
            0.851 0.791 0.808 0.856 0.900 0.914 0.974; ... % 2 Hz
            0.816 0.756 0.764 0.816 0.871 0.891 0.964; ... % 2.5 Hz
            0.773 0.721 0.730 0.776 0.841 0.866 0.953; ... % 3.15 Hz
            0.741 0.684 0.705 0.745 0.809 0.838 0.941; ... % 4 Hz
            0.726 0.628 0.678 0.736 0.780 0.812 0.929; ... % 5 Hz
            0.714 0.557 0.656 0.723 0.753 0.786 0.916; ... % 6.3 Hz
            0.670 0.520 0.623 0.678 0.728 0.765 0.904; ... % 8 Hz
            0.591 0.483 0.556 0.615 0.701 0.749 0.893; ... % 10 Hz
            0.554 0.446 0.523 0.614 0.685 0.737 0.884];    % 12.5 Hz
        
        % Settings
        AuditoryMasking = 2011;
        NoiseCorrection = 1;
        chans = 1;
        fs = nan; % needed for output structure
        fc = [125, 250, 500, 1000, 2000, 4000, 8000]; % not used directly
        FilterVersion = [];
        doplot = 1;
        
    else
        
        % WAVE FILE MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % in this mode a wav file is specified in the first input argument of
        % this function. A path can be concatenated to the file name.
        
        % read the wave file
        [data, fs]=audioread(IR);
        
        % If you have an older version of Matlab that does not support audioread,
        % then use the following instead:
        %[data, fs] = wavread(IR);
    end
else
    % DIRECT DATA ENTRY MODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % the first input argument is used for direct input of IR data
    data = IR;
    % required data from input structure
    if isstruct(IR)
        data = IR.audio;
        fs = IR.fs;
    end
    if nargin < 7
        fs = 48000; % default sampling rate
        disp('Audio sampling rate of 48 kHz has been assumed')
    end
end






% *************************************************************************
% LENGTH OF THE DATA AND NUMBER OF CHANNELS

if ~validate
    
    data=mean(data,3); % mix-down the third dimension if it exists
    
    [len,chans]=size(data);
    % Check length
    if len < fs*1.6
        disp('For STI calculation, the impulse response should be at least 1.6 s in duration.')
        disp('Zero padding has been applied, but the result should be interpreted with caution.')
        % zero padding is applied to the start of the IR
        data = [zeros(ceil(fs*1.6 - len),chans); data];
    end
    
    % Some extra zero-padding to avoid losing the temporal smearing of the octave
    % band filters in cases where the wave has substantial amplitude right at the
    % begining or end. 100 ms is added to each end.
    % Using the current filters, only zero padding at the tail is helpful,
    % but if zero phase filters are used then zero padding at the start
    % should also be helpful.
    data = [zeros(ceil(fs*0.1),chans); data; zeros(ceil(fs*0.1),chans)];
    [len,chans]=size(data); % update len
    
    % Note that zero padding may have a small effect on calculated mtfs,
    % especially when no noise or auditory effects are applied. For most
    % purposes, such effects are likely to be negligible.
    
end






% *************************************************************************
% DEFAULT SETTINGS

if ~validate
    % AARAE's linear phase filters by default
    if nargin < 8, FilterVersion = 2; end
    
    % do plot by default
    if nargin < 6, doplot = 1; end
    
    if ~dialogmode
        % do noise correction by default
        if nargin < 5, NoiseCorrection = 1; end
        
        % do auditory masking by default
        if nargin < 4, AuditoryMasking = 2011; end
    end
    
    if AuditoryMasking == 1, AuditoryMasking = 2011; end
    if AuditoryMasking == 2, AuditoryMasking = 2003; end
    
    % negligible noise by default
    if nargin < 3 && ~dialogmode
        Lnoise = repmat([-10 -10 -10 -10 -10 -10 -10], [chans,1]);
    else
        [row, col] = size(Lnoise);
        if col == 7 && row == 1, Lnoise = repmat(Lnoise, [chans,1]);
        elseif row == 7 && col ==1, Lnoise = repmat(Lnoise', [chans,1]);
        elseif ~(col == 7 && row == chans)
            warning('Lnoise is incorrectly dimensioned.')
            disp('Negligible background noise has been assumed.')
            Lnoise = repmat([-10 -10 -10 -10 -10 -10 -10], [chans,1]);
        end
    end
    
    
    % speech signal at 60 dBA by default
    % It is not advisable to use this default unless you have a good reason
    % to do so.
    if nargin < 2 && ~dialogmode
        Lsignal = 60 + [2.9, 2.9, -0.8, -6.8, -12.8, -18.8, -24.8];
        Lsignal = repmat(Lsignal, [chans,1]);
    else
        [row, col] = size(Lsignal);
        if col == 7 && row == 1, Lsignal = repmat(Lsignal, [chans,1]);
        elseif row == 7 && col ==1, Lsignal = repmat(Lsignal', [chans,1]);
        elseif ~(col == 7 && row == chans)
            warning('Lsignal is incorrectly dimensioned.')
            disp('Speech level of 60 dBA has been assumed.')
            Lsignal = 60 + [2.9, 2.9, -0.8, -6.8, -12.8, -18.8, -24.8];
            Lsignal = repmat(Lsignal, [chans,1]);
        end
    end
    
    if isfield(IR,'chanID')
        chanID = IR.chanID;
    else
        chanID = cellstr([repmat('Chan',size(data,2),1) num2str((1:size(data,2))')]);
    end
    
    % *************************************************************************
    % SET CONSTANTS AND PRE-ALLOCATE MATRICES
    
    
    % Nyquist frequency
    Nyquist=fs/2;
    
    % Time in seconds for each sample
    time=((1:len)-1)'./fs;
    
    
    % Define the octave band filter parameters for the basic filterbank
    bandnumber=21:3:39; % filter band numbers
    fc=10.^(bandnumber./10); % filter centre frequencies in Hz
    % Note that these filter settings are always used when auralization is done
    bandwidth = 1; % set to 0.5 instead of 1 if you want to try half-octave bandwidths
    f_low=fc./10^(0.15*bandwidth); % low cut-off frequency in Hz
    f_hi=fc.*10^(0.15*bandwidth); % high cut-off frequency in Hz
    
    halforder = 2; 
    
    b = zeros(halforder*2+1,length(fc)); % pre-allocate filter coefficients
    a = b; % pre-allocate filter coefficients
    
    % calculate filter coefficients
    % These filters are not class 1 or class 0 (passband is not flat
    % enough)
    for k = 1:length(fc)
        [b(:,k), a(:,k)]=butter(halforder, [f_low(k)/Nyquist f_hi(k)/Nyquist]);
    end
    
end

if doplot
    % constants used for plotting
    
    % screen position & size for figures
    bdwidth = 5;
    topbdwidth = 30;
    set(0,'Units','pixels')
    scnsize = get(0,'ScreenSize');
    pos1  = [bdwidth,...
        1/2*scnsize(4) + bdwidth,...
        scnsize(3)/2 - 2*bdwidth,...
        scnsize(4)/1.2 - (topbdwidth + bdwidth)];
    pos_step = 5;
    
    % rainbow colours
    colors = [255, 0, 0; ... % red
        255, 128, 0; ... % orange
        204, 204, 0; ... % dark yellow
        0, 204, 0; ... % mid green
        0, 204, 204; ... % dark cyan
        0, 0, 255; ... % blue
        127, 0, 255]; % violet
    colors = colors / 255; % rescale to 0-1 range
end

% List the modulation frequencies
% mf = [0.63 0.8 1 1.25 1.6 2 2.5 3.15 4 5 6.3 8 10 12.5]; % nominal frequencies
mf = 10.^((-2:11)/10); % exact frequencies

% pre-allocate outputs
M_STI = zeros(chans,1); % male STI
F_STI = zeros(chans,1); % female STI
STIPA = zeros(chans,1); % STIPA
MTI=zeros(chans,length(fc)); % modulation transfer indices
MTF = zeros(length(mf),length(fc),chans); % modulation transfer function

% pre-allocate outputs without noise correction nor auditory masking
M_STI0 = zeros(chans,1); % male STI
F_STI0 = zeros(chans,1); % female STI
STIPA0 = zeros(chans,1); % STIPA
STI_1985 = zeros(chans,1);
MTI0=zeros(chans,length(fc)); % modulation transfer indices
MTF0 = zeros(length(mf),length(fc),chans); % modulation transfer function






% *************************************************************************
% CHANNEL LOOP
% Values are calculated in the order of input channels

for ch = 1:chans
    
    
    %***************************************************************
    % Calculate Modulation Transfer Function (MTF) for 14 modulation
    % frequencies and 7 octave bands
    
    if ~validate
        mtfmethod = 1; % for testing purposes only. Set to ~=1 to try
        % FFT-based MTF analysis (see 'else...')
        
        % Modulation Transfer Function (MTF) calculations
        P_octave=zeros(length(data(:,ch)), length(length(fc)));
        MTF_ch=zeros(length(mf), length(fc));
        %MTFfilt_ch=zeros(length(mf), length(fc));
        if FilterVersion == 1
            P_octave = octavebandfilters(data(:,ch), fs);
        elseif FilterVersion == 2
            % use AARAE's linear phase filters  (recommended)
            % pseudo-Butterworth response
            orderin = 12; % in-band filter pseudo-order
            orderout = 12; % out-of-band filter pseudo-order
            if exist('FilterStrength','var')
                orderin = orderin * FilterStrength;
                orderout = orderout * FilterStrength;
            end
            P_octave = octbandfilter_viaFFT(data(:,ch),fs,...
                [125,250,500,1000,2000,4000,8000],[orderin,orderout],0,...
                1000,0,0,2);
            %octbandfilter_viaFFT(IN,fs,...
            %param,order,zeropad,...
            %minfftlenfactor,test,phasemode,base)
            if exist('FilterComp','var')
                if FilterComp == 1
                    deltatest = fftshift([1; zeros(size(data,1)-1,1)]);
                    deltatest_octave = octbandfilter_viaFFT(deltatest,fs,...
                [125,250,500,1000,2000,4000,8000],[orderin,orderout],0,...
                1000,0,0,2); 
                end
            end
        end
        for k=1:length(fc)
            if FilterVersion == 0
                % linear phase filter, but frequency selectivity does not
                % meet class 0 or 1
                P_octave(:,k)=filtfilt(b(:,k),a(:,k), data(:,ch)); 
            end
            if mtfmethod == 1
                for j=1:length(mf)
                    % number of whole number cycles to use for each modulation frequency
                    Fm_cycles = floor(len .* mf(j)./fs);
                    
                    % number of samples to use for each modulation frequency
                    Fm_len = floor(fs.*Fm_cycles ./ mf(j));
                    
                    % derive MTF using Schroeder's formula
                    
                    % direct implementation of calculation
                    MTF_num=abs(sum(P_octave(1:Fm_len,k).^2.*exp(-2i*pi*mf(j).*time(1:Fm_len))));
                    MTF_den=sum(P_octave(1:1:Fm_len,k).^2);
                    
                    MTF_ch(j,k)=MTF_num/MTF_den;
                    
                    if exist('FilterComp','var')
                        if FilterComp == 1
                            MTFfilt_num=abs(sum(deltatest_octave(1:Fm_len,k).^2.*exp(-2i*pi*mf(j).*time(1:Fm_len))));
                            MTFfilt_den=sum(deltatest_octave(1:1:Fm_len,k).^2);
                            %MTFfilt_ch(j,k)=MTFfilt_num/MTFfilt_den;
                            % Apply filter compensation here
                            MTF_ch(j,k) = MTF_ch(j,k) ./ (MTFfilt_num/MTFfilt_den);
                        end
                    end
                end
            else
                % FFT-derived modulation transfer function
                % This code is not optimised (it is not neccessary for it
                % to be in a for loop), and it was just added for
                % testing purposes.
                % This method is not recommended because there is unlikely
                % to be a good match between the modulation frequencies and
                % the fft components in the low frequency range. Using a
                % longer FFT solves this, but adds unnecessarily to
                % computational load.
                fftlen = len*1; % use a value > 1 (e.g. 10) to improve the
                % match between desired and actual frequencies. When this
                % number is large, the results are similar to the other
                % method.
                
                % magnitude spectrum of the squared envelope
                envspectrum = abs(fft(P_octave(:,k).^2,fftlen));
                
                % frequency vector
                f_env = fs * (0:(fftlen-1)) / fftlen;
                
                for j = 1:length(mf)
                    % find closest match to mf
                    ind = find(abs(mf(j)-f_env) == min(abs(mf(j)-f_env)),1,'first');
                    % modulation transfer ratio value
                    MTF_ch(j,k) = envspectrum(ind) ./ envspectrum(1);
                end
            end
        end
    end
    MTF0(:,:,ch) = MTF_ch;
    
    
    %***************************************************************
    % Calculate output parameters prior to noise correction and auditory
    % masking.
    [M_STI0(ch), F_STI0(ch), STIPA0(ch), MTI0(ch,:), STI_1985(ch)] = CalculateSTI(MTF_ch);
    
    
    
    %***************************************************************
    % Noise, masking and reception threshold adjustments
    
    % Intensity of noise
    if NoiseCorrection
        In = 10.^(Lnoise(ch,:)./10);
    else
        In = zeros(1,7);
    end
    
    % Intensity of speech
    Is = 10.^(Lsignal(ch,:)./10);
    
    % Intensity of auditory masking, and of auditory reception threshold
    switch AuditoryMasking
        case 2011
            [Iam, Irt] = AM(Lsignal(ch,:), Lnoise(ch,:), 2011);
        case 2003
            [Iam, Irt] = AM(Lsignal(ch,:), Lnoise(ch,:), 2003);
        case -2011 % Ed 4 auditory masking without reception threshold
            Iam = AM(Lsignal(ch,:), Lnoise(ch,:), 2011);
            Irt = zeros(1,7);
        case -2003 % Ed 4 auditory masking without reception threshold
            Iam = AM(Lsignal(ch,:), Lnoise(ch,:), 2003);
            Irt = zeros(1,7);
        case -1 % reception threshold without auditory masking
            [~, Irt] = AM(Lsignal(ch,:), Lnoise(ch,:), 2011);
            Iam = zeros(1,7);
        otherwise
            Iam = zeros(1,7);
            Irt = zeros(1,7);
    end
    
    % Apply ajustments to MTF
    for k = 1:7
        
        MTF_ch(:,k) = MTF_ch(:,k) * (Is(k)) ...
            / (Is(k) + Iam(k) + Irt(k) + In(k));
    end
    
    
    % Alternative way of writing the above
    %     if NoiseCorrection
    %         for k = 1:7
    %         MTF_ch(:,k)=MTF_ch(:,k) ./ (1+10.^((Lnoise(ch,k)-Lsignal(ch,k))/10));
    %         end
    %     end
    %
    %         % Intensity of auditory masking, and of auditory reception threshold
    %     switch AuditoryMasking
    %         case 2011
    %             [Iam, Irt, I] = AM(Lsignal(ch,:), Lnoise(ch,:), 2011);
    %         case 2003
    %             [Iam, Irt, I] = AM(Lsignal(ch,:), Lnoise(ch,:), 2003);
    %
    %     end
    %
    %     if AuditoryMasking>0
    %     for k = 1:7
    %         MTF_ch(:,k)=MTF_ch(:,k)*(I(k))...
    %                /(I(k)+Iam(k)+Irt(k));
    %     end
    %     end
    %
    
    
    
    
    % limiting the MTF values to 1 at the maximum.
    MTF_ch(MTF_ch>1) = 1;
    
    % add to the MTF output matrix
    MTF(:,:,ch) = MTF_ch;
    
    
    
    
    
    %***************************************************************
    % Calculate STI values, having applied the noise correction and masking
    % settings.
    
    [M_STI(ch), F_STI(ch), STIPA(ch), MTI(ch,:)] = CalculateSTI(MTF_ch);
    
    
    
    %***************************************************************
    % Plotting
    
    if doplot
        
        pos = pos1;
        pos(1) = pos(1) + (ch-1)*pos_step;
        % Create figure
        figure('name',['STI_IR Channel ', num2str(ch)],'position', pos);
        
        % plot of mtfs
        subplot(3,1,1)
        hold on
        
        % x axis
        set(gca,'XTickLabel',{'0.63', '0.8', '1', '1.25', '1.6', '2',...
            '2.5', '3.15','4','5','6.3','8','10','12.5'}, ...
            'XTick',1:14)
        xlim([1 14]);
        
        % plot the mtf for each octave band (DisplayName can be used if a
        % legend is added). Use rainbow colours as defined above.
        plot(MTF_ch(:,1),'Color',colors(1,:),'DisplayName','125 Hz');
        plot(MTF_ch(:,2),'Color',colors(2,:),'DisplayName','250 Hz');
        plot(MTF_ch(:,3),'Color',colors(3,:),'DisplayName','500 Hz');
        plot(MTF_ch(:,4),'Color',colors(4,:),'DisplayName','1 kHz');
        plot(MTF_ch(:,5),'Color',colors(5,:),'DisplayName','2 kHz');
        plot(MTF_ch(:,6),'Color',colors(6,:),'DisplayName','4 kHz');
        plot(MTF_ch(:,7),'Color',colors(7,:),'DisplayName','8 kHz');
        
        % Uncomment the following line to add a legend
        % legend('show');
        
        % y-axis limits
        ylim([0 1]);
        
        % Create xlabel
        xlabel('Modulation Frequency (Hz)');
        
        % Create ylabel
        ylabel('Modulation Transfer Ratio');
        
        % Create title
        if chans == 1
            title(['Male STI ', num2str(round(M_STI(ch)*100)/100), ...
                '; Female STI ', num2str(round(F_STI(ch)*100)/100), ...
                '  (AM=', num2str(AuditoryMasking), ...
                '; NC=', num2str(NoiseCorrection),')'])
        else
            title(['Channel ', num2str(ch), ...
                ': Male STI ', num2str(round(M_STI(ch)*100)/100), ...
                '; Female STI ', num2str(round(F_STI(ch)*100)/100), ...
                '  (AM=', num2str(AuditoryMasking), ...
                '; NC=', num2str(NoiseCorrection),')'])
        end
        hold off

        subplot(3,1,2)
        % Bar plot of modulation transfer indices
        
        % Plot the values
        bh=bar(MTI(ch,:));
        
        % The following code enables each bar to be coloured individually,
        % using the same colours as the mtf subplot (hence the bar plot
        % acts as a legend for the mtf plot).
        chil=get(bh,'children');
        cd=repmat(1:7,5,1);
        cd=[cd(:);nan];
        colormap(colors);
        set(chil,'facevertexcdata',cd);
        
        % x-axis
        set(gca,'XTickLabel',{'125', '250', '500', '1k', '2k', '4k', '8k'})
        xlabel('Octave Band Centre Frequency (Hz)')
        
        % y-axis
        ylabel('Modulation Transfer Index')
        ylim([0 1])
        
        % show the value of each bar
        for k = 1:7
            if MTI(ch,k) < 0.8
                % black number above the bar
                text(k-0.25,MTI(ch,k)+0.05, ...
                    num2str(round(MTI(ch,k)*1000)/1000),'Color','k')
            else
                % or white number below the top of the bar
                % (for high MTI values)
                text(k-0.25,MTI(ch,k)-0.05, ...
                    num2str(round(MTI(ch,k)*1000)/1000),'Color',[1 1 1])
            end
        end
        
        subplot(3,1,3)
        % Plot of octave band signal and noise levels
        hold on
        
        % x-axis
        set(gca,'XTickLabel',{'125', '250', '500', '1k', '2k', '4k', '8k'})
        xlabel('Octave Band Centre Frequency (Hz)')
        
        % y-axis
        ymax = 10*ceil(max(max([Lsignal;Lnoise]))/10);
        ylim([0 ymax])
        ylabel('SPL (dB)')
        
        % Sound pressure level of the speech signal
        plot(Lsignal(ch,:),'b','Marker','o','DisplayName','Signal')
        Levels(:,1,ch) = Lsignal(ch,:);
        % Sound pressure level of physical background noise
        plot(Lnoise(ch,:),'Color', [0.6 0.2 0], ...
            'Marker','o','DisplayName','Physical Noise')
        Levels(:,2,ch) = Lnoise(ch,:);
        if AuditoryMasking > 0
            
            % Sound pressure level of masking from the band below
            plot(10*log10(Iam),'Color',[0 0.5 0], ...
                'Marker','x','DisplayName','Masking')
            Levels(:,3,ch) = 10*log10(Iam);
            % Sound pressure level of the auditory reception threshold
            plot(10*log10(Irt),'k','LineStyle',':', ...
                'DisplayName','Threshold')
            Levels(:,4,ch) = 10*log10(Irt);
            % Sum of all sources of noise
            Inoise_sum = In + Irt + Iam;
            plot(10*log10(Inoise_sum), 'r', ...
                'LineWidth',1,'LineStyle','--', ...
                'Marker','+', 'DisplayName','Total Noise')
            Levels(:,5,ch) = 10*log10(Inoise_sum);
        end
        legend('show','Location','EastOutside');
        hold off
        
        if isstruct(IR)
            f = figure('Name',['STI_IR Channel ',num2str(ch),': STI, MTI & MTF'], ...
                'Position',[200 200 630 540]);
            %[left bottom width height]
            dat1 = MTF_ch;
            cnames1 = {'125', '250', '500', '1k', '2k', '4k', '8k'};
            rnames1 = {'0.63', '0.8', '1', '1.25', '1.6', '2',...
                '2.5', '3.15','4','5','6.3','8','10','12.5'};
            t1 =uitable('Data',dat1,'ColumnName',cnames1,'RowName',rnames1);
            %set(t,'ColumnWidth',{60});
            
            dat2 = [MTI(ch,:);Lsignal(ch,:);Lnoise(ch,:)];
            cnames2 = {'125', '250', '500', '1k', '2k', '4k', '8k'};
            rnames2 = {'MTI', 'Signal SPL (dB)', 'Noise SPL (dB)'};
            t2 =uitable('Data',dat2,'ColumnName',cnames2,'RowName',rnames2);
            %set(t,'ColumnWidth',{60});
            
            
            dat3 = [M_STI(ch);F_STI(ch);STIPA(ch);M_STI0(ch);F_STI0(ch);STIPA0(ch);STI_1985(ch)];
            cnames3 = {'value'};
            rnames3 = {'STI male','STI female','STIPA(IR)','M_STI without masking','F_STI0 without masking','STIPA without masking','STI 1985 version'};
            t3 =uitable('Data',dat3,'ColumnName',cnames3,'RowName',rnames3);
            %set(t,'ColumnWidth',{100});
            
            [~,tables] = disptables(f,[t3 t2 t1],{['Chan ' num2str(ch) ' - STI'],['Chan ' num2str(ch) ' - MTI'],['Chan ' num2str(ch) ' - MTF']});
            Verbose.tables = [Verbose.tables tables];
        end
    end % if doplot
    
    
end % channel loop
if isstruct(IR) && doplot >= 0
    mf = [0.63,0.8,1,1.25,1.6,2,2.5,3.15,4,5,6.3,8,10,12.5];
    doresultleaf(MTF,'Modulation transfer ratio',{'Modulation_frequency'},...
                 'Modulation_frequency', num2cell(mf),                                'Hz',          true,...
                 'Frequency',            num2cell([125,250,500,1000,2000,4000,8000]), 'Hz',          false,...
                 'Channel',              chanID,                                      'categorical', [],...
                 'name','Modulation_TF_STI');

    if AuditoryMasking > 0
        doresultleaf(Levels,'SPL [dB]',{'Frequency'},...
                     'Frequency', num2cell([125,250,500,1000,2000,4000,8000]),                   'Hz',          true,...
                     'Level',     {'Signal','Physical noise','Masking','Threshold','Total S+N'}, 'categorical', [],...
                     'Channel',   chanID,                                                        'categorical', [],...
                     'name','Band_SPL_STI');
    else
        doresultleaf(Levels,'SPL [dB]',{'Frequency'},...
                     'Frequency', num2cell([125,250,500,1000,2000,4000,8000]), 'Hz',          true,...
                     'Level',     {'Signal','Physical noise'},                 'categorical', [],...
                     'Channel',   chanID,                                      'categorical', [],...
                     'name','Band_SPL_STI');
    end
end

%***************************************************************
% In the case of single channel analysis, remove the singleton dimension

if chans == 1
    M_STI = squeeze(M_STI);
    F_STI = squeeze(F_STI);
    MTF = squeeze(MTF);
end







%***************************************************************
% Create verbose output structure

% outputs for the MTF calculated as per input arguments
Verbose.M_STI = M_STI;
Verbose.F_STI = F_STI;
Verbose.STIPA = STIPA;
Verbose.MTI = MTI;
Verbose.MTF = MTF;

% outputs replicating the signal and noise inputs
Verbose.Lsignal = Lsignal;
Verbose.Lnoise = Lnoise;

% outputs for the MTF without any treatment for noise or masking
Verbose.M_STI0 = M_STI0;
Verbose.F_STI0 = F_STI0;
Verbose.STIPA0 = STIPA0;
Verbose.STI_1985 = STI_1985;
Verbose.MTI0 = MTI0;
Verbose.MTF0 = MTF0;

% output text reporting the function settings

if ischar(IR)
    inputstring = IR;
else
    inputstring = 'Direct IR input.';
end

Verbose.txt = [datestr(now), '  ', inputstring, '  fs:', num2str(fs), ...
    '  Auditory masking:', num2str(AuditoryMasking), ...
    '  Noise correction:', num2str(NoiseCorrection), ...
    '  Filterbank:', num2str(FilterVersion)];

Verbose.funcallback.name = 'STI_IR.m';
Verbose.funcallback.inarg = {fs,Lsignal,Lnoise,AuditoryMasking,NoiseCorrection,doplot,FilterVersion,FilterStrength,FilterComp};

varargout{1} = M_STI;
varargout{2} = F_STI;


%***************************************************************
% Simple auralization for 1 and 2 chan analysis (only available in dialogmode)

if doAuralization
    % Use a dialog box to select a wav file.
    % The variable speechrec is used to store the file name.
    [speechrec, pathname] = uigetfile({'*.wav'}, ...
        'Open a wav file containing a dry speech recording for auralization');
    
    % Open the wav file, and acquire IR data and sampling rate.
    [speechwave, fs2] = audioread([pathname,speechrec]);
    
    % If you have an older version of Matlab that does not support audioread,
    % then use the following instead:
    %[speechwave, fs2] = wavread([pathname,speechrec]);
    
    % Alternatively, a dry speech recording could be hard-coded to load
    % here, thereby avoiding the need for a dialog box.
    
    speechwave = speechwave(:,1); % use channel 1 only
    
    % match sampling rate
    if ~(fs2 == fs)
        speechwave = resample(speechwave,fs,fs2);
    end
    
    lenspeech = length(speechwave);
    
    % an even number of samples is required for the frequency domain
    % operations that follow
    outputlength = 2*ceil((len + lenspeech - 1)/2);
    
    % number of IR channels to auralize (maximum of 2)
    if chans > 2
        AuralizeChans = 2;
    else
        AuralizeChans = chans;
    end
    
    % Convolve speech with IR by multiplying spectra
    auralization = ifft(fft(repmat(speechwave,[1,AuralizeChans]), outputlength) ...
        .* fft(data(:,1:AuralizeChans), outputlength));
    
    
    
    % Adjust speech spectrum in each octave band **************************
    
    % Filter convolved speech into octave bands
    S_octave = zeros(outputlength,AuralizeChans,7); % pre-allocate
    for k=1:length(fc)
        S_octave(:,:,k)=filter(b(:,k),a(:,k), auralization); % filter
    end
    
    % Get the unadjusted rms octave band level
    S_level0 = (10*log10(mean(S_octave.^2)));
    
    % Add 3 dB, based on Annex J.4 of the standard
    % (we are assuming that Lsignal has already been determined based on
    % Annex J or similar procedures/assumptions).
    S_level0 = S_level0 + 3;
    
    % Adjust the level of the speech auralization in each octave band
    S_octave = S_octave .* ...
        repmat(10.^((reshape(Lsignal(1:AuralizeChans,:),1,AuralizeChans,7) ...
        - S_level0)./20), [outputlength,1,1]);
    
    
    
    
    % Generate noise ******************************************************
    
    % First generate pink noise using a frequency domain approach:
    
    % Magnitude slope function (for half spectrum, not including DC and
    % Nyquist)
    magslope = ((1:outputlength/2-1)./(outputlength/4)).^(-0.5)';
    
    % low freq noise cutoff (80 Hz)
    lowcut = floor(80 * outputlength / fs);
    magslope(1:lowcut) = 0; % zero the values below cutoff
    
    % generate noise in the frequency domain, by random phase
    noisyslope = repmat(magslope,1,AuralizeChans) ...
        .* exp(1i*2*pi.*rand(outputlength/2-1,AuralizeChans));
    clear magslope;
    
    % Transform to time domain
    PinkNoise = ifft([zeros(1,AuralizeChans);noisyslope; ...
        zeros(1,AuralizeChans);flipud(conj(noisyslope))]);
    clear noisyslope;
    
    % Filter the noise into octave bands
    N_octave = zeros(outputlength,AuralizeChans,7); % pre-allocate
    for k=1:length(fc)
        N_octave(:,:,k)=filter(b(:,k),a(:,k), PinkNoise); % filter
    end
    
    % Get the unadjusted rms octave band level
    N_level0 = 10*log10(mean(N_octave.^2));
    
    % Adjust the level of the noise auralization in each octave band
    N_octave = N_octave .* ...
        repmat(10.^((reshape(Lnoise(1:AuralizeChans,:),1,AuralizeChans,7) ...
        - N_level0)./20), [outputlength,1,1]);
    
    
    
    % Final steps in the auralization *************************************
    
    % Combine the noise and speech (summing the octave band filtered
    % signals)
    auralization = sum(S_octave,3) + sum(N_octave,3);
    
    % Normalize the auralization.
    % Note that this code could be improved later to allow for calibrated
    % auralization. Without calibrating the auralization, the auditory
    % threshold and masking effects are not re-created properly.
    auralization = auralization ./ max(max(abs(auralization)));
    
    % Play the auralization.
    sound(auralization, fs)
    
    % Loop for replaying, saving and finishing
    choice = 'x'; % create a string
    
    % loop until the user presses the 'Done' button
    while ~strcmp(choice,'Done')
        choice = questdlg('What next?', ...
            'Auralization ...', ...
            'Play again', 'Save audio', 'Done','Done');
        switch choice
            case 'Play again'
                sound(auralization, fs)
            case 'Save audio'
                [filename, pathname] = uiputfile({'*.wav'},'Save as');
                audiowrite([pathname,filename],auralization,fs);
        end % switch
    end % while
    
end % if doAuralization


end % eof








%***************************************************************
%***************************************************************
% Functions for auditory masking,
% STI calculation and filterbuilder octave band filters
%***************************************************************
%***************************************************************



%***************************************************************
function [Iam, Irt, I]=AM(signal, noise, version)
% Auditory masking function

% 2011 version of the standard is default
if nargin < 3, version = 2011; end

% signal level plus noise level
level_combined=10*log10(10.^(signal./10)+10.^(noise./10));

% Intensity of level_combined in octave bands
I=10.^(level_combined./10);

% preallocate
amf=zeros(1,7);     % Auditory Masking Factor
amfdB=zeros(1,7); % Auditory Masking Factor in dB
Iam=zeros(1,7);   % Masking Intensity

switch version
    case 2011
        % auditory masking from 2011 version of the standard
        for k=2:7
            if level_combined(k-1)< 63
                amfdB(k)=0.5*level_combined(k-1)-65;
            elseif level_combined(k-1) < 67;
                amfdB(k)=1.8*level_combined(k-1)-146.9;
            elseif level_combined(k-1) < 100;
                amfdB(k)=0.5*level_combined(k-1)-59.8;
            else
                amfdB(k)=-10;
            end
            amf(k)=10^(amfdB(k)/10);
            Iam(k)=I(k-1)*amf(k);
        end
        
    case 2003
        % auditory masking from 2003 version of the standard
        for k=2:7
            if round(level_combined(k-1))>=46 && round(level_combined(k-1)) <=55
                amf(k)=0.0001;
            elseif round(level_combined(k-1))>=56 && round(level_combined(k-1)) <=65;
                amf(k)=0.000316;
            elseif round(level_combined(k-1)) >=66 && round(level_combined(k-1)) <=75;
                amf(k)=0.003162;
            elseif round(level_combined(k-1)) >=76 && round(level_combined(k-1)) <=85;
                amf(k)=0.01;
            elseif round(level_combined(k-1)) >=86 && round(level_combined(k-1)) <=95;
                amf(k)=0.031622;
            elseif round(level_combined(k-1)) >95;
                amf(k)=0.1;
            end
            
            Iam(k)=I(k-1)*amf(k);
        end
end

% Absolute Speech Reception Threshold
ART=[46 27 12 6.5 7.5 8 12];

% Intensity of threshold
Irt=10.^(ART./10);
end % eof








%***************************************************************
function [M_STI, F_STI, STIPA, MTI, STI_1985] = CalculateSTI(MTF)
% calculate STI values and modulation transfer indices
% from a modulation transfer function matrix

% Effecive signal to noise ratio (SNR)
SNReff=10*log10(MTF./(1-MTF));

% limit values to -15 <= effSNR <= 15 dB
SNReff(SNReff>15)=15;
SNReff(SNReff<-15)=-15;

% Houtgast 1985
SNmean = mean(SNReff,1);
wk = [0.13,0.14,0.11,0.12,0.19,0.17,0.14];
SNmean = sum(SNmean.*wk);
STI_1985 = (SNmean+15)./30;

% Calculate Transmission Index (TI) - 2003 & 2011
% and averaged Modulation Tranasfer Index (MTI)
TI=(SNReff+15)./30;
MTI=mean(TI,1);


% STI Male
alpha=[0.085 0.127 0.230 0.233 0.309 0.224 0.173];
beta=[0.085 0.078 0.065 0.011 0.047 0.095];
MTI_alpha=alpha.*MTI;
MTI_beta=beta.*sqrt(MTI(1:length(beta)).*MTI(2:length(beta)+1));
M_STI=sum(MTI_alpha)-sum(MTI_beta);


% STI Female
alpha=[0 0.117 0.223 0.216 0.328 0.250 0.194];
beta=[0 0.099 0.066 0.062 0.025 0.076];
MTI_alpha=alpha.*MTI;
MTI_beta=beta.*sqrt(MTI(1:length(beta)).*MTI(2:length(beta)+1));
F_STI=sum(MTI_alpha)-sum(MTI_beta);


% STIPA(IR) (2011 version)
alpha=[0.085 0.127 0.230 0.233 0.309 0.224 0.173];
beta=[0.085 0.078 0.065 0.011 0.047 0.095];
MTI_STIPA = [mean([TI(5,1); TI(12,1)]), ...
    mean([TI(3,2); TI(10,2)]), ...
    mean([TI(1,3); TI(8,3)]), ...
    mean([TI(6,4); TI(13,4)]), ...
    mean([TI(4,5); TI(11,5)]), ...
    mean([TI(2,6); TI(9,6)]), ...
    mean([TI(7,7); TI(14,7)])];
MTI_alpha=alpha.*MTI_STIPA;
MTI_beta=beta.*sqrt(MTI_STIPA(1:length(beta)) ...
    .*MTI_STIPA(2:length(beta)+1));
STIPA=sum(MTI_alpha)-sum(MTI_beta);

% % STI 1985 version
% wk = [0.13 0.14 0.11 0.12 0.19 0.17 0.14];
% STI_1985 = (sum(wk.* MTI) + 15)/30;
% disp(num2str(STI_1985))
end % eof





function y = octavebandfilters(x, fs)
% Octave band filters generated from Matlab's filterbuilder are applied to
% the input wave, x.
% (This may not be compatible with some versions of Matlab.)

y = zeros(length(x), 7); % preallocate y
B  = 1;       % Bands per octave
N  = 6;       % Order



warning('off') % avoids warnings about F0 centre frequency



F0 = 10^(21/10);  % Center frequency
h125 = fdesign.octave(B, 'Class 0', 'N,F0', N, F0, fs);
Hd125 = design(h125, 'butter', ...
    'SOSScaleNorm', 'Linf');



F0 = 10^(24/10);
h250 = fdesign.octave(B, 'Class 0', 'N,F0', N, F0, fs);
Hd250 = design(h250, 'butter', ...
    'SOSScaleNorm', 'Linf');



F0 = 10^(27/10);
h500 = fdesign.octave(B, 'Class 0', 'N,F0', N, F0, fs);
Hd500 = design(h500, 'butter', ...
    'SOSScaleNorm', 'Linf');



F0 = 1000;
h1000 = fdesign.octave(B, 'Class 0', 'N,F0', N, F0, fs);
Hd1000 = design(h1000, 'butter', ...
    'SOSScaleNorm', 'Linf');



F0 = 10^(33/10);
h2000 = fdesign.octave(B, 'Class 0', 'N,F0', N, F0, fs);
Hd2000 = design(h2000, 'butter', ...
    'SOSScaleNorm', 'Linf');



F0 = 10^(36/10);
h4000 = fdesign.octave(B, 'Class 0', 'N,F0', N, F0, fs);
Hd4000 = design(h4000, 'butter', ...
    'SOSScaleNorm', 'Linf');



F0 = 10^(39/10);
h8000 = fdesign.octave(B, 'Class 0', 'N,F0', N, F0, fs);
Hd8000 = design(h8000, 'butter', ...
    'SOSScaleNorm', 'Linf');




y(:,1) = filter(Hd125, x);
y(:,2) = filter(Hd250, x);
y(:,3) = filter(Hd500, x);
y(:,4) = filter(Hd1000, x);
y(:,5) = filter(Hd2000, x);
y(:,6) = filter(Hd4000, x);
y(:,7) = filter(Hd8000, x);

warning('on')

end % eof

