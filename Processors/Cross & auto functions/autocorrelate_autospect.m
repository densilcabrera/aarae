function OUT = autocorrelate_autospect(in)
% Performs autocorrelation by multiplication in the frequency domain
% (autospectrum method). The output is normalized.
%
% Code by Densil Cabrera
% version 1.0 (11 October 2013)

len = length(in.audio); % length of the IR

% fft length must be double the signal length to avoid circular autocorrelation
fftlen = len * 2;

% derive spectrum
spectrum = fft(in.audio,fftlen);

% multiply spectrum with its conjugate, and return to time domain
out = ifft(spectrum .* conj(spectrum));

% normalize
out = out ./ repmat(out(1,:,:,:,:,:),[fftlen,1,1,1,1,1]);

if isstruct(in)
    OUT.audio = out;
    OUT.funcallback.name = 'autocorrelate_autospect.m';
    OUT.funcallback.inarg = {};
else
    OUT = out;
end

