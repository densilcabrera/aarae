function OUT = Aweight(IN,fs)
% A-weighting using Matlab's Audio Toolbox
if nargin < 2
    if isstruct(IN)
        audio = IN.audio;
        fs = IN.fs;
    else
        audio = IN;
        fs = inputdlg({'Sampling frequency [samples/s]'},...
            'Fs',1,{'48000'});
        fs = str2num(char(fs));
    end
else
    audio = IN;
end
if ~isempty(audio) && ~isempty(fs)
    weightType = 'A-weighting';
    weightFilt = weightingFilter(weightType,fs);
    audio = weightFilt(audio);
    if isstruct(IN)
        OUT = IN;
        OUT.audio = audio;
        OUT.funcallback.name = 'Aweight.m';
        OUT.funcallback.inarg = {};
    else
        OUT = processed;
    end
end