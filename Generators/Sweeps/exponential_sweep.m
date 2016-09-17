function OUT = exponential_sweep(dur,start_freq,end_freq,fs,reverse,rcos_ms)% generates an exponentially swept
% This function generates an exponential sweep (also known as logarithmic
% sweep), which is one of the most commonly used test signals for measuring
% impulse responses.


% signal S, starting at start_freq Hz and ending at end_freq Hz,
% for duration = dur seconds long, and an
% amplitude of ampl = 0.5, with a raised cosine window applied for rcos_ms = 15 ms.
% Sinv is the inverse of S.
if nargin == 0
    param = inputdlg({'Duration [s]';...
                       'Start frequency [Hz]';...
                       'End frequency [Hz]';...
                       'Sampling Frequency [samples/s]';...
                       'Ascending [0] or descending [1] sweep';...
                       'Fade-in and Fade-out duration [ms]'},...
                       'Sine sweep input parameters',1,{'10';'20';'20000';'48000';'0';'15'});
    param = str2num(char(param));
    if length(param) < 6, param = []; end
    if ~isempty(param)
        dur = param(1);
        start_freq = param(2);
        end_freq = param(3);
        fs = param(4);
        reverse = param(5);
        rcos_ms = param(6);
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

    OUT.audio = S';
    OUT.audio2 = Sinv';
    OUT.fs = fs;
    OUT.tag = ['Sine sweep exp' num2str(dur)];
    OUT.properties.dur = dur;
    OUT.properties.sig_len = sig_len;
    OUT.properties.freq = [start_freq, end_freq];
    OUT.properties.reverse = reverse;
    OUT.properties.rcos_ms = rcos_ms;
    OUT.properties.scale_inv = scale_inv;
    OUT.properties.IRscalingfactor = IRscalingfactor; % used by convolveaudiowithaudio2.m
    OUT.funcallback.name = 'exponential_sweep.m';
    OUT.funcallback.inarg = {dur,start_freq,end_freq,fs,reverse,rcos_ms};
else
    OUT = [];
end