function OUT = pure_tone(duration, frequency, fs)
% Generates a single channel pure tone at the specified frequency
%
% The tone is a cosine wave (a sine wave is also generated in audio2)
if nargin == 0
    param = inputdlg({'Duration [s]';...
                       'Center frequency [Hz]';...
                       'Sampling Frequency [samples/s]';...
                       'Phase offset [radians]'},...
                       'Sine sweep input parameters',1,...
                       {'10';'1000';'48000';'0'});
    param = str2num(char(param));
    if length(param) < 4, param = []; end
    if ~isempty(param)
        duration = param(1);
        frequency = param(2);
        fs = param(3);
        phi = param(4);
    end
end
if ~isempty(param) || nargin ~= 0
    t = linspace(0,duration,fs*duration);
    y = cos(2*pi*frequency.*t'+phi);
    y2 = sin(2*pi*frequency.*t'+phi);
    tag = [num2str(frequency) 'Hz tone'];
    
    OUT.audio = y;
    OUT.audio2 = y2;
    OUT.fs = fs;
    OUT.tag = tag;
    OUT.funcallback.name = 'pure_tone.m';
    OUT.funcallback.inarg = {duration, frequency, fs};
else
    OUT = [];
end

end % End of function