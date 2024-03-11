function OUT = Cweight(IN,fs)
% C weighting using Matlab's Audio Toolbox
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
    weightType = 'C-weighting';
    weightFilt = weightingFilter(weightType,fs);
    audio = weightFilt(audio);
    if isstruct(IN)
        OUT = IN;
        OUT.audio = audio;
        OUT.funcallback.name = 'Cweight.m';
        OUT.funcallback.inarg = {};
    else
        OUT = audio;
    end
end