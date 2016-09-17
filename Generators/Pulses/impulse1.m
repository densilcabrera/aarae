function OUT = impulse1(duration_pre,duration_post,amplitude, fs)
% This function generates a Kronecker delta function with a specified
% duration of silence before and/or after. This provides a signal that can
% be used for the direct measurement of impulse response, although it is
% highly susceptible to contamination from noise.
%
% Using the AARAE generator dialog box, before running this function, you
% can generate a pulse train, which can be useful for:
%  * Synchronous averaging (to remove background noise contamination); or
%  * Examining time variance in low noise systems.
%
% As well as the audio wave, this function outputs its inverse filter
% (which is trivial - i.e. a simple time reversal of the delta function) so
% that AARAE's convolve audio with audio2 (*) button can be used to perform
% synchronous averaging or IR stacking (for time variance analysis) from a
% recording of a pulse train. Problems may occur if the gap between
% impulses when generated in the pulse train is too short for the system's
% latency (however a gap of 1 s should usually be safe, and delay
% compensation in AARAE's Syscal may help).
%
% Code by Densil Cabrera
% version 2 (1 August 2014)

if nargin == 0
    param = inputdlg({'Duration of silence before impulse [s]';...
        'Duration of silence after impulse [s]';...
        'Impulse amplitude value';...
        'Sampling frequency [samples/s]'}, ...
        'Impulse input parameters',1,{'0';'0.1';'1';'48000'});
    param = str2num(char(param));
    if length(param) < 4, param = []; end
    if ~isempty(param)
        duration_pre = param(1);
        duration_post = param(2);
        amplitude = param(3);
        fs = param(4);
    end
else
    param = [];
end
if ~isempty(param) || nargin ~= 0
    samples_pre = round(duration_pre * fs);
    samples_post = round(duration_post * fs);
    
    if samples_pre > 0
        wave_pre = zeros(samples_pre,1);
    else
        wave_pre = [];
    end
    
    if samples_post > 0
        wave_post = zeros(samples_post,1);
    else
        wave_post = [];
    end
    
    y = [wave_pre; amplitude; wave_post];
    tag = ['Impulse_' num2str(samples_pre) '_' num2str(samples_post) '_' num2str(amplitude)];
    
    OUT.audio = y;
    OUT.audio2 = flipud(y);
    OUT.fs = fs;
    OUT.tag = tag;
    OUT.funcallback.name = 'impulse1.m';
    OUT.funcallback.inarg = {duration_pre,duration_post,amplitude, fs};
else
    OUT = [];
end

end % End of function