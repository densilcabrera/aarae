function OUT = Cmessageweight(IN,fs)
% C message weighting
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
        content = load([cd '/Processors/Filters/' num2str(fs) 'Hz/Cmessage-WeightingFilter.mat']);
        filterbank = content.filterbank;
        processed = filter(filterbank,1,audio);
        
    else
        % Matlab's filterbuilder
        WT = 'Cmessage';  % Weighting type
        
        
        h = fdesign.audioweighting('WT', WT, fs);
        
        Hd = design(h, 'bell41009', 'SOSScaleNorm', 'Linf');
        processed = filter(Hd,audio);
    end
    if isstruct(IN)
        OUT = IN;
        OUT.audio = processed;
        OUT.funcallback.name = 'Cmessageweight.m';
        OUT.funcallback.inarg = {};
    else
        OUT = processed;
    end
end