function OUT = ITUR4684weight(IN,fs)
% ITU-R weighting for audio frequency noise measurement
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
        content = load([cd '/Processors/Filters/' num2str(fs) 'Hz/ITUR4684-WeightingFilter.mat']);
        filterbank = content.filterbank;
        processed = filter(filterbank,1,audio);
        
    else
        % Matlab's filterbuilder
        WT = 'ITUR4684';  % Weighting type
        
        
        h = fdesign.audioweighting('WT', WT, fs);
        
        Hd = design(h, 'iirlpnorm', 'SOSScaleNorm', 'Linf');
        processed = filter(Hd,audio);
    end
    if isstruct(IN)
        OUT = IN;
        OUT.audio = processed;
        OUT.funcallback.name = 'ITUR4684weight.m';
        OUT.funcallback.inarg = {};
    else
        OUT = processed;
    end
end