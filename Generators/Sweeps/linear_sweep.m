% Generates a linear sweep and its inverse for IR measurement
%
function OUT = linear_sweep(dur,start_freq,end_freq,reverse,fs)


if nargin == 0
    param = inputdlg({'Duration [s]';...
                       'Start frequency [Hz]';...
                       'End frequency [Hz]';...
                       'Sampling Frequency [samples/s]';...
                       'Ascending [0] or descending [1] sweep'},...
                       'Sine sweep input parameters',1,{'10';'0';'24000';'48000';'0'});
    param = str2num(char(param));
    if length(param) < 5, param = []; end
    if ~isempty(param)
        dur = param(1);
        start_freq = param(2);
        end_freq = param(3);
        fs = param(4);
        reverse = param(5);
    end   
else
    param = [];
end
if ~isempty(param) || nargin ~=0
    if (exist('fs') ~= 1)
       fs = 48000;
    end
    
    if end_freq > fs/2, end_freq = fs/2; end
    
    t = 0:1/fs:(dur-1/fs); %time in seconds
    %S = chirp(t,start_freq,(dur-1/fs),end_freq)';
    t1 = dur-1/fs;
    p = 1;
    phi = 0;
    rcos_ms = 15;
    beta   = (end_freq-start_freq).*(t1.^(-p));
    S = sin(2*pi * ( beta./(1+p).*(t.^(1+p)) + start_freq.*t + phi/360));
    rcos_len = round(length(S)*((rcos_ms*1e-3)/dur));
    sig_len = length(S);
    rcoswin = hann(2*rcos_len).';
    S = [S(1:rcos_len).*rcoswin(1:rcos_len),S(rcos_len+1:sig_len-rcos_len),S(sig_len-rcos_len+1:sig_len).*rcoswin((rcos_len+1):(rcos_len*2))]';
    if reverse == 1
        Sinv = S;
        S = flipud(S);
    else
        Sinv = flipud(S);
    end
    
    % IR scaling factor (this could be done much more simply)
    fftS = fft(S);
    Sinvfft = fft(Sinv);
    mid_freq = (start_freq + end_freq)/2;
    index = round(mid_freq/(fs/sig_len));
    const1 = abs(conj(fftS(index))/(abs(fftS(index))^2));
    const2 = abs(Sinvfft(index));
    IRscalingfactor = const1/const2;
    

    OUT.audio = S;
    OUT.audio2 = Sinv;
    OUT.fs = fs;
    OUT.tag = ['Sine sweep linear' num2str(dur)];
    OUT.properties.dur = dur;
    OUT.properties.freq = [start_freq,end_freq];
    OUT.properties.reverse = reverse;
    OUT.properties.IRscalingfactor = IRscalingfactor; % used by convolveaudiowithaudio2.m
    OUT.funcallback.name = 'linear_sweep.m';
    OUT.funcallback.inarg = {dur,start_freq,end_freq,reverse,fs};
else
    OUT = [];
end