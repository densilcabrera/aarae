function OUT = Bweight(IN,fs)
% B-weighting filter (calling code from Mike Brookes' Voicebox)

% Note that stdspectrum (from Voicebox) is in aarae's Utilities folder.
processed = [];
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
    %if isdir([cd '/Processors/Filters/' num2str(fs) 'Hz'])
    if false % bypass ths code
        content = load([cd '/Processors/Filters/' num2str(fs) 'Hz/B-WeightingFilter.mat']);
        filterbank = content.filterbank;
        processed = filter(filterbank,1,audio);
        
    else
        % use stdspectrum.m from voicebox
        [b,a] = stdspectrum(3,'z',fs);
        processed = filter(b,a,audio);
    end
    if isstruct(IN)
        OUT = IN;
        OUT.audio = processed;
        OUT.funcallback.name = 'Bweight.m';
        OUT.funcallback.inarg = {};
    else
        OUT = processed;
    end
end