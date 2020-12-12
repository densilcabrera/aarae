function [OUT, varargout] = Silence(fs, duration, chans, bands)
% This function generates silence (i.e. a waveform of zeros).
%
% Note that multi-band signals are not compatible with playback (use 1 band
% for signals that can be played).


if nargin == 0 % If the function is called within the AARAE environment it
               % won't have any input arguments, this is when the inputdlg
               % function becomes useful.

    param = inputdlg({'Sampling rate (Hz)';... 
                      'Duration (s)'; ...
                      'Number of channels';...
                      'Number of bands'},...
                      'Settings',... % This is the dialog window title.
                      [1 60],... 
                      {'48000';'1';'1';'1'}); % Default values
    param = str2num(char(param)); 

    if length(param) < 4, param = []; end 
    if ~isempty(param) 
        fs = param(1);
        duration = param(2);
        chans = param(3);
        bands = param(4);
    end
else
    param = [];
end


if ~isempty(param) || nargin ~= 0
    

    
    audio = zeros(round(duration*fs),chans,bands);
    
    
    if nargin == 0
        OUT.audio = audio; 
        OUT.fs = fs;
        OUT.tag = 'Silence'; 
        OUT.properties.generator = 'Silence';
        OUT.funcallback.name = 'Silence.m';
        OUT.funcallback.inarg = {fs, duration, chans, bands};
    end
    
   
    if nargin ~= 0
        OUT = audio;
        varargout{1} = fs;
    end
else
    % AARAE requires that in case that the user doesn't input enough
    % arguments to generate audio to output an empty variable.
    OUT = [];
end

end % End of function