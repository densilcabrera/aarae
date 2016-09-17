function [OUT, varargout] = InverseHilbert(IN,fs)
% This function performs an inverse Hilbert transform (only useful for
% complex input data)


    
    
if isstruct(IN) 
    audio = IN.audio; % Extract the audio data
    fs = IN.fs;       % Extract the sampling frequency of the audio data
else
    audio = IN;
    if nargin < 2
        fs = inputdlg({'Sampling frequency [samples/s]'},...
                           'Fs',1,{'48000'});
        fs = str2num(fs);
    end
end

if ~isempty(audio)
   % inverse of the Hilbert transform
    audio = abs(audio).* cos(angle(audio));
    if isstruct(IN)
        OUT = IN; 
        OUT.audio = audio;
        OUT.funcallback.name = 'InverseHilbert.m';
        OUT.funcallback.inarg = {fs};
    else
        OUT = audio;
    end
    varargout{1} = fs;
else
    OUT = [];
end

