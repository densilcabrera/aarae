function OUT = audio_shuffle(in,fs,windowlength,depthfactor, outduration, reverse, playaudio)
% Fragments and shuffles an audio waveform to create a new audio texture.
% The process is a bit like granular or concatenative synthesis.
% This can be useful as a sonification tool.
%
% Loosely based on the processes described in:
% S. Ferguson and D.Cabrera (2009) "Auditory spectral summarisation for 
% audio signals with musical applications," 10th International Society for 
% Music Information Retrieval Conference, Kobe, Japan, 567-572.
%
% Code by Densil Cabrera
% version 1.0 (12 October 2013)
%
% INPUT ARGUMENTS
% in is the audio data (can be multi-channel or single channel). It can be
% a structure with the required fields .audio and .fs - in which case the
% other parameters are set via a dialog box.
%
% fs is audio sampling rate in Hz
%
% windowlength is the mean window length (in samples)of each audio window 
% to be shuffled. In dialog box mode, windowlength is controlled indirectly
% by entering the window time in seconds.
%
% depthfactor is the number of windows x windowlength ./ number of samples
% in the original data. Greater depthfactor results in a more constant
% sound texture.
%
% outduration is the output wave duration in seconds
%
% reverse can be set or clear (0 | 1): if it is set then approximately half
% of the windows will be time reversed; if it is clear then all of the
% windows will be forward in time.

if isstruct(in)
    in = choose_from_higher_dimensions(in,2,1); 
    % required field of input structure
    data = in.audio; % audio waveform
    fs = in.fs; % audio sampling rate
else
    % if 'in' is not a structure, then it is an audio waveform
    data = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end    
end
if nargin < 7, playaudio = 0; end
if nargin < 6, reverse = 1; end
if nargin < 5, outduration = 1; end
if nargin < 4, depthfactor = 4; end
if nargin < 3,
    windowlength = .125;
    % dialog box for settings
     prompt = {'Window duration (s):', ...
        'Depth factor:', ...
        'Output duration (s):', ...
        'Include time reverse (0 | 1):', ...
        'Play sound (0 | 1)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'0.125','4','1','1','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(answer)
        OUT = [];
        return
    else
        windowlength = str2num(answer{1,1})*fs;
        depthfactor = str2num(answer{2,1});
        outduration = str2num(answer{3,1});
        reverse = str2num(answer{4,1});
        playaudio = str2num(answer{5,1});
    end 
end

if ~isempty(data) && ~isempty(fs) && ~isempty(windowlength) && ~isempty(depthfactor) && ~isempty(outduration) && ~isempty(reverse) && ~isempty(playaudio)
    data = squeeze(data(:,:,1)); % remove 3rd dimension if it exists
    [len, chan] = size(data);

    % output length in samples
    outlen = outduration * fs;

    % pre-allocate output
    y = zeros(outlen, chan);

    numberofwindows = ceil(depthfactor * outlen ./ windowlength);
    disp(['Number of windows in audio shuffle: ', num2str(numberofwindows)])

    maxwinlength = ceil(windowlength * 2.^0.5);
    if maxwinlength > 0.2 * len
        % keep maximim window length to no more than 20% of the audio data
        % length
        maxwinlength = 0.2 * len;
        disp(['Maximum window length has been reduced to ', num2str(maxwinlength/fs)])
    end
    % minwinlength = floor(windowlength ./ 2.^0.5);



    for n = 1:numberofwindows
        % determine data range and get data for the given window
        specificwinlen = round(windowlength * 2.^(rand(1,1) - 0.5));
        numberofavailablewindowcentres = len - specificwinlen - 2;   
        windowcentre = round(rand(1,1) * numberofavailablewindowcentres + specificwinlen./2 +1);
        startindex = windowcentre-floor(specificwinlen./2);
        endindex = startindex + specificwinlen-1;
        wf = window(@hann, specificwinlen); 
        datawindow = data(startindex:endindex,:) .* repmat(wf,[1 chan]);

        if reverse && rand(1)> 0.5
            datawindow = flipud(datawindow);
        end

        % write data to output
        numberofavailablewindowcentres = outlen - specificwinlen - 2; 
        windowcentre = round(rand(1,1) * numberofavailablewindowcentres + specificwinlen./2 +1);
        startindex = windowcentre-floor(specificwinlen./2);
        endindex = startindex + specificwinlen-1;
        y(startindex:endindex,:) = y(startindex:endindex,:) + datawindow;

    end
    %normalize
    y = y ./max(max(abs(y)));

    if isstruct(in)
        OUT = in;
        OUT.audio = y;
        OUT.funcallback.name = 'audio_shuffle.m';
        OUT.funcallback.inarg = {fs,windowlength,depthfactor,outduration,reverse,playaudio};
    else
        OUT = y;
    end
    
    if playaudio
        sound(y,fs)

            % Loop for replaying, saving and finishing
        choice = 'x'; % create a string

        % loop until the user presses the 'Done' button
        while ~strcmp(choice,'Done')
            choice = questdlg('What next?', ...
                'Audio Play ...', ...
                'Play again', 'Save audio', 'Done','Done');
            switch choice
                case 'Play again'
                    sound(y, fs)
                case 'Save audio'
                    [filename, pathname] = uiputfile({'*.wav'},'Save as');
                    audiowrite([pathname,filename],y,fs);
            end % switch
        end % while
    end
else
    OUT = [];
end