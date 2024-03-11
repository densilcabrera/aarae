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
    for d3 = 1:size(audio,3)
        for d4 = 1:size(audio,4)
            for d5 = 1:size(audio,5)
                for d6 = 1:size(audio,6)
                    audio(:,:,d3,d4,d5,d6) = weightFilt(audio(:,:,d3,d4,d5,d6));
                end
            end
        end
    end
    if isstruct(IN)
        OUT = IN;
        OUT.audio = audio;
        OUT.funcallback.name = 'Cweight.m';
        OUT.funcallback.inarg = {};
    else
        OUT = audio;
    end
end