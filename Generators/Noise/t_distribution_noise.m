function OUT = t_distribution_noise(degreesoffreedom,duration, chans, fs)
% This function generates random numbers from the Student's t distribution.
% A very large degreesoffreedom value yields a distribution approaching
% Gaussian. When degreesoffreedom == 1, the function generates Cauchy
% noise, which is likely to have a very large range of values.

if nargin == 0
    param = inputdlg({'Degrees of freedom (integer >=1)';...
        'Duration of the wave [s]';...
        'Number of channels';...
        'Sampling frequency [samples/s]';...
        'Normalize the waveform? [0 | 1]'}, ...
        'Input parameters',1,{'1';'1';'1';'48000';'1'});
    param = str2num(char(param));
    if length(param) < 5, param = []; end
    if ~isempty(param)
        degreesoffreedom = round(param(1));
        duration = param(2);
        chans = round(param(3));
        fs = param(4);
        normalize = param(5);
    end
else
    param = [];
end

if degreesoffreedom <1, degreesoffreedom = 1; end
if chans <1, chans = 1; end

if ~isempty(param) || nargin ~= 0
    y = trnd(degreesoffreedom,round(duration*fs),chans);
    if normalize == 1
        y = y ./max(max(abs(y)));
    end

    tag = ['t_noise' num2str(duration)];
    
    OUT.audio = y;
    %OUT.audio2 = flipud(y); % inverse that assumes flat spectrum (poor assumption)
    OUT.audio2 = ifft(1./fft(y)); % alternative 'inverse'
    OUT.fs = fs;
    OUT.tag = tag;
    OUT.funcallback.name = 't_distribution_noise.m';
    OUT.funcallback.inarg = {degreesoffreedom,duration,chans,fs,normalize};
else
    OUT = [];
end

end % End of function