function OUT = pure_tone_complex(duration, frequency, fs)
% Generates a single channel pure tone at the specified frequency, with a
% real and imaginary component
%
% The resulting wave is complex-valued.
%
% Note that this function does NOT generate a harmonic series tone (the
% term 'complex' is not used in that sense).

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
    y = cos(2*pi*frequency.*t'+phi)+1i*sin(2*pi*frequency.*t'+phi);
    tag = [num2str(frequency) 'Hz tone'];
    
    OUT.audio = y;
    OUT.fs = fs;
    OUT.tag = tag;
    OUT.funcallback.name = 'pure_tone_complex.m';
    OUT.funcallback.inarg = {duration, frequency, fs};
else
    OUT = [];
end

end % End of function