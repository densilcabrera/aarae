function OUT = ExponentialEnvelope(in, fs, exponent, normalize, starttime)
% This function applies an exponential envelope to the audio.
%
% Use a positive exponent for exponential growth.
% Use a negative exponent for exponential decay.
%
% For equivalent reverberation times, the following can be used:
% RT    | exponent
% 0.5 s | -13.8
% 1 s   | -6.9
% 2 s   | -3.45
%
% In other words, the exponent is equal to -6ln(10) / (2RT)
%
%
% Code by Densil Cabrera
% version 1.00 (4 November 2013)

if nargin < 4, normalize = 1; end
if nargin < 3
    prompt = {'Exponential growth or decay constant', ...
        'Normalize', 'Start time (s)'};
    dlg_title = 'Settings';
    num_lines = 1;
    def = {'0','1','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        OUT = [];
        return
    else
        exponent = str2double(answer{1,1});
        normalize = abs(str2double(answer{2,1}));
        starttime = str2double(answer{3,1});
    end
end
if isstruct(in)
    audio = in.audio;
    fs = in.fs;
else
    audio = in;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
end

if ~isempty(audio) && ~isempty(fs) && ~isempty(exponent) && ~isempty(normalize)
    [len, chans, bands,dim4,dim5,dim6] = size(audio);
    
    startind = abs(round(starttime*fs)) + 1;
    if startind > len
        startind = 1;
    end
    t = (0:(len-startind))./fs; % time in s
    envelope = exp(exponent.*t)';
    
    audio(startind:end,:,:,:,:,:) = audio(startind:end,:,:,:,:,:)...
        .* repmat(envelope, [1,chans,bands,dim4,dim5,dim6]);
    if normalize
        audio = audio ./ max(max(max(max(max(max(abs(audio)))))));
    end
    
    if isstruct(in)
        OUT.audio = audio;
        OUT.funcallback.name = 'ExponentialEnvelope.m';
        OUT.funcallback.inarg = {fs,exponent,normalize,starttime};
    else
        OUT = audio;
    end

    RT = -3*log(10)/exponent;
    disp(['Equivalent reverberation time of exponent is ', num2str(RT), ' s.'])
else
    OUT = [];
end