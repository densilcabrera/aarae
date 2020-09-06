function OUT = complementarysweeps(dur,gap,start_freq,end_freq,fs,reverse,rcos_ms)
% This function generates a pair of sweeps, where the second sweep is the
% same as the first multiplied by -1. This pair of test signals may be
% useful for suppressing tonal noise (e.g. electrical hum) by making the
% second sweep start a whole number of periods (of the tonal noise) after 
% the first sweep.
%
% For example, if we wish to suppress a 50 Hz mains power hum, we can
% separate the sweeps by a multiple of 1/50 = 0.02 s. Actually, simply
% separating them by a whole number of seconds fulfills that.
%
% Currently this generator only generates exponential sweeps.



if nargin == 0
    param = inputdlg({'Duration [s]';...
                       'Silence between sweeps [s]';...
                       'Start frequency [Hz]';...
                       'End frequency [Hz]';...
                       'Sampling Frequency [samples/s]';...
                       'Ascending [0] or descending [1] sweep';...
                       'Fade-in and Fade-out duration [ms]'},...
                       'Sine sweep input parameters',1,{'10';'1';'20';'20000';'48000';'0';'15'});
    param = str2num(char(param));
    if length(param) < 7, param = []; end
    if ~isempty(param)
        dur = abs(param(1));
        gap = abs(param(2));
        start_freq = abs(param(3));
        end_freq = abs(param(4));
        fs = abs(param(5));
        reverse = param(6);
        rcos_ms = param(7);
    end
else
    param = [];
end
if ~isempty(param) || nargin ~=0
    if ~exist('fs','var')
       fs = 48000;
    end
    SI = 1/fs;
    ampl = 0.5;
    %rcos_ms = 15;
    scale_inv = 1;
    if 2*round((rcos_ms*1e-3))>dur, rcos_ms = 15; end

    w1 = 2*pi*start_freq; w2 = 2*pi*end_freq;
    K = (dur*w1)/(log(w2/w1));
    L = log(w2/w1)/dur;
    t = [0:round(dur/SI)-1]*SI;
    phi = K*(exp(t*L) - 1);
    freq = K*L*exp(t*L);
    %freqaxis = freq/(2*pi);
    amp_env = 10.^((log10(0.5))*log2(freq/freq(1)));
    S = ampl*sin(phi);
    rcos_len = round(length(S)*((rcos_ms*1e-3)/dur));
    sig_len = length(S);
    rcoswin = hann(2*rcos_len).';
    S = [S(1:rcos_len).*rcoswin(1:rcos_len),S(rcos_len+1:sig_len-rcos_len),S(sig_len-rcos_len+1:sig_len).*rcoswin((rcos_len+1):(rcos_len*2))];
    Sinv = fliplr(S).*amp_env;

    % correction for allpass delay
    Sinvfft = fft(Sinv);
    Sinvfft = Sinvfft.*exp(1i*2*pi*(0:(sig_len-1))*(sig_len-1)/sig_len);
    Sinv = real(ifft(Sinvfft));

    if scale_inv == 1
       fftS = fft(S);
       mid_freq = (start_freq + end_freq)/2;
       index = round(mid_freq/(fs/sig_len));
       const1 = abs(conj(fftS(index))/(abs(fftS(index))^2));
       const2 = abs(Sinvfft(index));
       IRscalingfactor = const1/const2;
       % Sinv = Sinv * IRscalingfactor; % scaling factor is applied in
       % convolveaudiowithaudio2
    else
        IRscalingfactor = 1;
    end
    
    if reverse
        S = fliplr(S);
        Sinv = fliplr(Sinv);
    end

    OUT.audio = [S';zeros(round(fs*gap),1);-S';zeros(round(fs*gap),1)];
    OUT.audio2 = Sinv';
    OUT.fs = fs;
    OUT.tag = ['ComplSweeps' num2str(dur)];
    OUT.properties.dur = dur;
    OUT.properties.complementarysweepsgap = gap;
    OUT.properties.sig_len = sig_len;
    OUT.properties.freq = [start_freq, end_freq];
    OUT.properties.reverse = reverse;
    OUT.properties.rcos_ms = rcos_ms;
    OUT.properties.scale_inv = scale_inv;
    OUT.properties.IRscalingfactor = IRscalingfactor; % used by convolveaudiowithaudio2.m
    OUT.funcallback.name = 'complementarysweeps.m';
    OUT.funcallback.inarg = {dur,gap,start_freq,end_freq,fs,reverse,rcos_ms};
else
    OUT = [];
end
